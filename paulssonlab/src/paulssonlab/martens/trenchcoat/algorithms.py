from params import parse_segmentation_channel_names
from properties import write_properties_to_table, make_cell_type

# TODO: abstract out the grouping / iteration kind (i.e. could segment on pairs of channels etc.)
# But for now, just make a rudimentary version.
def run_segmentation_analysis(h5file, name, fov, frame, seg_params, channels, seg_algorithm, out_dir_masks, out_dir_tables, regions, cell_dict):
    # Get the channel names for segmentation
    seg_channels = parse_segmentation_channel_names(seg_params)
    
    # Create directory structure to store the masks
    dir_path = "{}/File_{}/FOV_{}/".format(out_dir_masks, name, fov)
    pathlib.Path(dir_path).mkdir(parents=True, exist_ok=True)
    h5file_masks = tables.create_file("{}/File_{}/FOV_{}/Frame_{}.h5".format(out_dir_masks, name, fov, frame), mode="w")
    
    # Create directory structure to store the tables
    dir_path = "{}/File_{}/FOV_{}/".format(out_dir_tables, name, fov)
    pathlib.Path(dir_path).mkdir(parents=True, exist_ok=True)
    h5file_tables = tables.create_file("{}/File_{}/FOV_{}/Frame_{}.h5".format(out_dir_tables, name, fov, frame), mode="w")
    
    # Stack of image regions, or whole image?
    if regions:
        loop_func = make_ch_to_img_stack
    else:
        loop_func = make_ch_to_img
    
    # For writing tables
    Cell = make_cell_type(cell_dict)
    table = h5file_tables.create_table("/tables", Cell)
    row = table.row
    
    # Iterate all Z levels & all seg channels
    for z in h5file.iter_nodes("/Images/File_{}/FOV_{}/Frame_{}".format(name, fov, frame))
        z_level = z._v_name
        
        # Store images, or image stacks, by channel
        # Required for writing all intensity values at once, in the same table row.
        ch_to_img = loop_func(h5file, z, channels)
        
        # Masks & tables
        for ch in seg_channels:
            # Calculate the mask(s)
            masks = seg_algorithm(ch_to_img[ch], seg_params[ch])
            
            # Write the masks & tables
            for region_number, mask in enumerate(masks):
                # Masks
                h5file_masks.create_array("/{}/{}".format(z_level, ch), "region_{}".format(region_number), obj=mask, createparents=True)
                
                # Compute the properties & write to disk
                properties = skimage.measure.regionprops(mask)
                write_properties_to_table(h5file_tables, name, fov, frame, region_number, properties, params, row, ch_to_img]
    
    # Done!
    table.flush()
    h5file_masks.close()
    h5file_tables.close()

# Make a stack?
# Not None: analyze regions within the image
# Load stack of image slices with pre-determined dimensions.
# Ranges is a numpy array with (min_row, min_col, max_row, max_col) for each region.
# NOTE: Pass in the node reference, not the whole image, into this sub-routine.
# Then, when using the slice notation below, will take advantage of
# the chunking so as to not load the entire image.
def make_ch_to_img_stack(h5file, z, channels):
    ch_to_img = {}
    
    for ch in channels:
        image_node = h5file.get_node(z, ch)
        image = load_stack(image_node, regions)
        
        x_dimension = regions[0, 2] - regions[0, 0]
        y_dimension = regions[0, 3] - regions[0, 1]
        stack_size = regions.shape[0]
        
        stack = numpy.empty(shape=(stack_size, y_dimension, x_dimension), dtype=numpy.uint16)
        
        # min_row, min_col, max_row, max_col
        for i, r in enumerate(regions):
            stack[i] = image_node[r[0] : r[2], r[1] : r[3]]
        
        ch_to_img[c] = image
    
    return ch_to_img

# None: analyze the whole image
def make_ch_to_img(h5file, z, channels):
    ch_to_img = {}
    
    for ch in channels:
        image_node = h5file.get_node(z, ch)
        image = image_node.read()
        # Add extra z-axis to make it compatible with generic code which iterates over z
        image = image[numpy.newaxis]
        ch_to_img[ch] = image
        
    return ch_to_img

# Very basic thresholding!
def run_thresholding(stack, params):
    return (stack > params["threshold"]).astype(numpy.uint16)    

# Thresholding which also uses Niblack (mean, std dev) & watershedding
# Stack must have 3 dimensions
def run_niblack_segmentation(stack, params):
    if len(stack.shape) != 3:
        pass # ERROR!
    
    # Store a stack of labeled segmentation masks
    result = numpy.empty(stack.shape, dtype=numpy.uint16)
    
    # For each image in the stack
    for z, stack_elem in enumerate(stack):
        # Threshold is calculated using Otsu method, which samples the distribution of pixel intensities, and then
        # scaled by the otsu multiplier, which might differ depending on empirical experience with a given channel or dataset.
        threshold = skimage.filters.threshold_otsu(stack_elem) * params["otsu_multiplier"]
        
        # If it's a low value, then it's probably just noise...
        if threshold < params["garbage_otsu_value"]:
            result[z] = numpy.zeros(shape=(stack.shape[1], stack.shape[2]))
        else:
            # Set pixels < threshold to zero, but preserve values > threshold.
            # Helps with the mean + stdev*k calculations in Niblack.
            # It is "bimodal" because the values are either zero, or they are the original value.
            # The "multiply by scaling_factor" trick helps to emphasize the value of the pixels which are kept.
            # This enhances the std. dev. at the edges of cells.
            bimodal = (stack_elem < threshold) * 0 + \
                (stack_elem >= threshold) * stack_elem.astype(numpy.uint32) * params["scaling_factor"]

            # Apply Niblack method
            niblack = skimage.filters.threshold_niblack(bimodal, window_size=params["niblack_w"], k=params["niblack_k"])

            # Only keep pixels which are < than the mean*stdev*k (in rectangular region, size w).
            mask = (stack_elem < niblack)

            # Closing helps to fill in internal holes inside the cell
            mask = skimage.morphology.binary_closing(mask)

            # Remove small objects
            skimage.morphology.remove_small_objects(mask, min_size=params["min_size"], connectivity=2, in_place=True)

            # Label the regions
            markers = skimage.measure.label(mask)

            # Make a basin (invert max, min pixels)
            # NOTE it would be nice if the watershed algo. allowed for a max priority queue.
            # Would save us a step here.
            image = stack_elem.max() - stack_elem

            result[z] = skimage.morphology.watershed(image=image, markers=markers, mask=mask)
        
    return result

