# Input path to H5 file
# Return all channel names, for all files, as a dictionary of ndarrays
def get_all_channel_names(h5file_path):
    h5file = tables.open_file(h5file_path, mode="r")
    
    file_to_channels = {}
    
    for node in h5file.iter_nodes("/Metadata"):
        file_to_channels[node._v_name] = h5file.get_node(node, "channels")[:]
    
    h5file.close()
    
    return file_to_channels

# Return the channels just for a specific ND2 file (already opened)
def get_channel_names(h5file):
    return h5file.get_node(h5file.root, "channels")[:]

# Input path to parameters YAML file, return dictionary of parameters
def read_params_file(file_path):
    yaml = YAML(typ="safe", pure=True)
    return yaml.load(Path(file_path))

# Write dictionary of parameters to a YAML file
def write_params(params, file_path):
    yaml = YAML(typ="safe", pure=True)
    yaml.dump(params, Path(file_path))

# Input params dict, return list of fluorescence channels
def parse_fluorescence_channel_names(params):
    fluor_channels = []
    for channel in params.keys():
        if (params[channel]['fluorescent'] == True):
            fluor_channels.append(key)
    
    return fluor_channels

# Input params dict, return list of segmentation channels
def parse_segmentation_channel_names(params):
    segmentation_channels = []
    for channel in params.keys():
        if (params[channel]['is_segmentation'] == True):
            segmentation_channels.append(key)
    
    return segmentation_channels

# Input a node in the HDF5 metadata section, return a dict. with the following metadata
def get_metadata(n):
    metadata = { 'channels'       : n.channels.read()
                 'fields_of_view' : n.fields_of_view.read()
                 'frames'         : n.frames.read()
                 'width'          : n.width.read()
                 'height'         : n.height.read()
                 'pixel_microns'  : n.pixel_microns.read()
               }
    
    return metadata
