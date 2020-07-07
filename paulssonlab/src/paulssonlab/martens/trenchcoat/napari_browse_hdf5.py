import tables
import numpy
import napari
from dask import delayed
import dask.array

from params import read_params_file

"""    
NOTE for now, assumes image dtype is float32

TODO clean up the metadata functions, and their names, then stick them in the separate metadata Python file
Might also need to rework the similar functions used for nd2 browsing?
"""

def get_largest_extents_hdf5(h5file, metadata_key):
    """
    TODO: for now, assume that the metadata types are integers,
    but in the future, could explicitly check, and set the min, max values accordingly.
    
    TODO: Also, support non list type objects.
    """
    from sys import maxsize
    
    smallest_value = maxsize
    largest_value = 0
    
    # only iterate the File_{} nodes
    metadata_node = h5file.get_node("/Metadata")()
    for node in h5file.iter_nodes(metadata_node):
        values = h5file.get_node(node, metadata_key).read()
        
        if values[0] < smallest_value:
            smallest_value = values[0]
        if values[-1] > largest_value:
            largest_value = values[-1]
    
    return (smallest_value, largest_value)

def metadata_attributes_equal(h5file, attribute):
    """
    Input an h5file with metadata, and an attribute of interest
    Returns None if not all nd2 files (in the h5file) have identical attributes of this type in their metadata
    Returns the attribute if they are identical
    """
    # only iterate the File_{} nodes
    metadata_node = h5file.get_node("/Metadata")()
    iter_nodes = h5file.iter_nodes(metadata_node) 
    zeroth_attribute = h5file.get_node(next(iter_nodes), attribute).read()
    
    for node in iter_nodes:
        next_attribute = h5file.get_node(node, attribute).read()
        if next_attribute != zeroth_attribute:
            return None
    
    return zeroth_attribute

def metadata_array_equal(h5file, attribute):
    """
    Input an h5file with metadata, and an attribute of interest
    Returns None if not all nd2 files (in the h5file) have identical attributes of this type in their metadata
    Returns the attribute if they are identical
    TODO Update this desription. Is for comparing arrays!
    """
    # only iterate the File_{} nodes
    metadata_node = h5file.get_node("/Metadata")()
    iter_nodes = h5file.iter_nodes(metadata_node) 
    zeroth_attribute = h5file.get_node(next(iter_nodes), attribute).read()
    
    for node in iter_nodes:
        next_attribute = h5file.get_node(node, attribute).read()
        if not numpy.array_equal(next_attribute, zeroth_attribute):
            return None
    
    return zeroth_attribute

def load_img(node, channel, height, width, dtype):
    """
    Attempt to load an image, and if it cannot be found, return a zero'ed array instead.
    """
    try:
        return node._f_get_child(channel).read()
    
    except:
        # FIXME width,height?
        return dask.array.from_delayed(numpy.zeros(shape=(height, width), dtype=dtype))

def main_hdf5_browser_function(images_file, masks_file, regions_file, napari_settings_file):
    """
    
    """
    # Init napari
    with napari.gui_qt():
        viewer = napari.Viewer()
        
        # Load the image layer params for napari
        # TODO params for masks layers?
        layer_params = read_params_file(napari_settings_file)
        
        # Define a function for lazily loading a single image from the h5 file
        lazy_load_img = delayed(load_img)
        
        # File with images & metadata
        h5file = tables.open_file(images_file, "r")
        
        # Get the largest extents across all the nd2 files within the h5 file,
        # so that the dask array is made with the proper min, max dimensions for each dimension.
        attributes = ['fields_of_view', 'frames', 'z_levels']
        extents = {}
        for a in attributes:
            (smallest, largest) = get_largest_extents_hdf5(h5file, a)
            if (smallest == largest):
                extents[a] = [smallest]
            else:
                extents[a] = [i for i in range(smallest, largest+1)]
        
        # Get the height & width, while checking that they are identical across all nd2 files
        height = metadata_attributes_equal(h5file, 'height')
        width = metadata_attributes_equal(h5file, 'width')
        
        # Define the channels
        channels = metadata_array_equal(h5file, 'channels')
        
        # Iterate the H5 file images & lazily load them into a dask array
        parent_node = h5file.get_node("/Images")
        file_nodes = [x._v_name for x in h5file.list_nodes(parent_node)]
        
        for c in channels:
        c = c.decode("utf-8")
        
        # File nodes
        # FIXME what order to iterate? do they need to be sorted?
        file_array = []
        for file in file_nodes:
            # FOV nodes
            fov_array = []
            for fov in extents['fields_of_view']:
                # Frame nodes
                frame_array = []
                for frame in extents['frames']:
                    # Z nodes
                    z_array = []
                    for z in extents['z_levels']:
                        path = "/Images/{}/FOV_{}/Frame_{}/Z_{}".format(file, fov, frame, z)
                        node = h5file.get_node(path)
                        arr = dask.array.from_delayed(lazy_load_img(node, c, height, width, numpy.float32),
                                                        shape=(height, width), # FIXME or width,height?
                                                        dtype=numpy.float32)
                        z_array.append(arr)
                    
                    z_array_to_stack = dask.array.stack(z_array)
                    frame_array.append(z_array_to_stack)
                
                frame_array_to_stack = dask.array.stack(frame_array)
                fov_array.append(frame_array_to_stack)
            
            fov_array_to_stack = dask.array.stack(fov_array)
            file_array.append(fov_array_to_stack)
        
        megastack = dask.array.stack(file_array)
        viewer_func(megastack, **params[c])
        
        # Done loading image layers.
        
        ## ### Masks ###
        
        
        # NOTE it might be that the looping below is redundant with the looping above --
        # does it make sense to merge the loops?
        # 
        # Would require many repeated checks on whether or not masks & trenches were specified,
        # versus repeating the same loops again.
        # 
        # I think separating them is cleaner. Can make it less ugly by sticking things into subroutines?
        
        # Load segmentation masks?
        # because each region gets its own masks, they need to be merged!
        # this means that it is necessary to specify which regions were used to do so during the segmentation step.
        # TODO figure out how to require regions? (or, the same region will just be over-written, and the user will figure it out?)
        # As with the trenches, these regions the masks within these regions need to be combined into a single, large image.
        # However, here we an additional choice: whether or not each of the regions gets unique labels,
        # or if the labels are re-used across regions.
        # I like the idea of re-using them, because then the colors will correspond to the same numbering schemes.
        if masks_file:
            h5file_masks = tables.open_file(masks_file, "r")
            
            # If there were regions, then the masks need to be agglomerated into a larger image
            if regions_file:
                # FIXME how to determine the shape? height,width or width,height?
                zeroed_canvas = numpy.zeros(shape=(height, width), dtype=numpy.uint16)
                
                h5file_regions = tables.open_file(regions_file, "r")
                
                # FIXME best way to iterate the segmentation channels?
                # could either pass in the seg params file again,
                # or try to infer it from the data.
                # I prefer the idea of inferring it, since there it is. But do we just pick
                # the first available node?
                
                # 1. open the masks file
                # 2. define a recursive iterator which specifically looks for arrays
                # 3. call next() on the iterator to get the first available array
                # 4. examine the array's _v_name to determine what its seg channel was
                # ... but that only gives us the first seg channel.
                # Way to iterate at a specific level, and then get just the children of the first node found?
                # FIXME wait, don't we skip writing the mask if it's all zeroes? That means we can't query
                # the file for all the seg. channels, unless we scan the *entire* file & get all unique names.
                # can we add a metadata section to the masks file, where we embed the YAML file used as seg settings?
                # this way, we always know what the settings were for a given dataset!
                
                # 
                # 
                # 
                # FIXME how to determine the number of region sets? should each region set get its own labels layer?
                # that would mean sc x region_sets total number of layers.
                # 
                # NOTE that if the regions within a set overlap, then they will overwrite each other
                # when it comes time to merging them together in the same canvas.
                for sc in seg_channels:
                    
                    # File nodes
                    # FIXME what order to iterate? do they need to be sorted?
                    # Needs to be the same order as the images were loaded.
                    file_array = []
                    for file in file_nodes:
                        # FOV nodes
                        fov_array = []
                        for fov in extents['fields_of_view']:
                            # Frame nodes
                            frame_array = []
                            for frame in extents['frames']:
                                # Regions are always shared across Z, so it's OK
                                # to read them in now.
                                path = "/{}/FOV_{}/Frame_{}"
                                regions_coords = h5file_regions.get_node(path).read()
                                
                                # Z nodes
                                z_array = []
                                for z in extents['z_levels']:
                                    
                                    # TODO now, need to loop all the available regions,
                                    # and then copy over the corresponding data from the masks into
                                    # a single, "zeroed out" canvas
                                    
                                    path = "/{}/FOV_{}/Frame_{}/Z_{}/{}/".format(file, fov, frame, z, )
                                    masks_iter = h5file_masks.iter_nodes(path)
                                    
                                    # FIXME use zipped iterators?
                                    for i, (coords, masks) in enumerate(zip(regions_coords, masks_iter)):
                                        
                                        
                                        
                                    
                                    
                                    
                                    
                                    
                                    path = "/{}/FOV_{}/Frame_{}/Z_{}/{}/0/region_0".format(file, fov, frame, z, sc)
                                    node = h5file_regions.get_node(path)
                                    arr = dask.array.from_delayed(delayed(node.read()),
                                                                  # If there are no regions, then it must fill the entire
                                                                  # dimensions of the image.
                                                                  shape=(height, width), # FIXME or width,height?
                                                                  dtype=numpy.uint16)
                                    
                                    z_array.append(arr)
                                
                                z_array_to_stack = dask.array.stack(z_array)
                                frame_array.append(z_array_to_stack)
                            
                            frame_array_to_stack = dask.array.stack(frame_array)
                            fov_array.append(frame_array_to_stack)
                        
                        fov_array_to_stack = dask.array.stack(fov_array)
                        file_array.append(fov_array_to_stack)
                    
                    # TODO what params to pass in? use a yaml file?
                    megastack = dask.array.stack(file_array)
                    viewer.add_labels(megastack)
                
                
            
            
            
            ### For reference:
            ### h5file_masks.create_carray("/{}/{}/{}".format(z_level, sc, region_set_number),
            
            # Just load the first detected region information
            # FIXME best way to iterate the segmentation channels?
            # could either pass in the seg params file again,
            # or try to infer it from the data.
            # I prefer the idea of inferring it, since there it is. But do we just pick
            # the first available node?
            
            # No regions file
            else:
                
                for sc in seg_channels:
                    
                    # File nodes
                    # FIXME what order to iterate? do they need to be sorted?
                    # Needs to be the same order as the images were loaded.
                    file_array = []
                    for file in file_nodes:
                        # FOV nodes
                        fov_array = []
                        for fov in extents['fields_of_view']:
                            # Frame nodes
                            frame_array = []
                            for frame in extents['frames']:
                                # Z nodes
                                z_array = []
                                for z in extents['z_levels']:
                                    # If no regions, then default to the zeroth region set number and region
                                    path = "/{}/FOV_{}/Frame_{}/Z_{}/{}/0/region_0".format(file, fov, frame, z, sc)
                                    node = h5file.get_node(path)
                                    arr = dask.array.from_delayed(delayed(node.read()),
                                                                  # If there are no regions, then it must fill the entire
                                                                  # dimensions of the image.
                                                                  shape=(height, width), # FIXME or width,height?
                                                                  dtype=numpy.uint16)
                                    
                                    z_array.append(arr)
                                
                                z_array_to_stack = dask.array.stack(z_array)
                                frame_array.append(z_array_to_stack)
                            
                            frame_array_to_stack = dask.array.stack(frame_array)
                            fov_array.append(frame_array_to_stack)
                        
                        fov_array_to_stack = dask.array.stack(fov_array)
                        file_array.append(fov_array_to_stack)
                    
                    # TODO what params to pass in? use a yaml file?
                    megastack = dask.array.stack(file_array)
                    viewer.add_labels(megastack)
            
        
        # Load trench masks?
        # idea here is that the trench masks are stored as rectangular regions, so they need to be
        # converted into an "image" with labeled regions.
        # simplest way is to init. an empty numpy array, and then iterate each of the regions & fill in
        # the corresponding areas with the value corresponding to the trench (region) number.
        if regions_file:
            pass
    
    # Done!
    h5file.close()
