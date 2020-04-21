import tables

# Once a given trench has a segmentation mask, then write the properties of the masked region to the table.
# NOTE future versions of skimage will allow spitting out a properties object all at once, rather than lazily calculating them one at a time.
def write_properties_to_table(filename, fov, frame_number, region_number, properties, params, row, sc_enum, channel_to_img):
    # NOTE skimage is changing their coordinates -- do we still want to transpose??? I think so...
    coords = properties.coords.T
    
    # Background corrections
    for ch, img in channel_to_img:
        row['total_intensity_{}'.format(ch)] = \
            subtract_background_from_coords(coords, img, params[ch]['background'])
    
    row['info_filename']      = filename
    row['info_fov']           = fov
    row['info_frame']         = frame_number
    row['info_seg_channel']   = sc_enum
    row['info_region_number'] = region_number
    row['info_label']         = properties.label
    
    row['geometry_Area']        = properties.area
    row['geometry_Orientation'] = properties.orientation
    row['geometry_Perimeter']   = properties.perimeter
    
    # FIXME very small numbers --> math domain error.
    # Shouldn't the skimage library handle this issue and just return NaN?
    # The values might drop below zero due to floating-point error,
    # and then a square root of a neg number causes the error.
    # Maybe we can wrap each of these in try clauses, and write NaN if they throw exceptions.
    row['axis_length_major'] = properties.major_axis_length
    row['axis_length_minor'] = properties.minor_axis_length
    
    row['centroid_row'] = properties.centroid[0]
    row['centroid_col'] = properties.centroid[1]
    
    row['bounding_box_min_row'] = properties.bbox[0]
    row['bounding_box_min_col'] = properties.bbox[1]
    row['bounding_box_max_row'] = properties.bbox[2]
    row['bounding_box_max_col'] = properties.bbox[3]
    
    # Append the properties information to the table
    row.append()    

# Given a dictionary, convert into a PyTables Table cell type
def make_cell_type(column_types):
    # Yeah, this syntax is really goofy
    return type('Cell', (tables.IsDescription,), column_types)

# Define the column types for the PyTable: this stores segmented cell information
# Returns a dictionary, which can then be used to make the type.
# NOTE: the type cannot be pickled or sent through processes, hence the 2 steps.
def make_cell_dict(channels, seg_channels, file_names):
    # Define the PyTables column types using a dictionary
    column_types = {
        "info_fov"             : tables.UInt16Col(),
        "info_frame"           : tables.UInt16Col(),
        
        # e.g. trench number. The "region" within the whole image.
        "info_region_number"   : tables.UInt16Col(),
        
        # The labeled, connected component within the region
        "info_label"           : tables.UInt16Col(),
        
        "geometry_Area"        : tables.UInt16Col(),
        "geometry_Orientation" : tables.Float32Col(),
        "geometry_Perimeter"   : tables.Float32Col(),
        
        # See note below about math domain errors.
        "axis_length_major"    : tables.Float32Col(),
        "axis_length_minor"    : tables.Float32Col(),
        
        "centroid_row"         : tables.UInt16Col(),
        "centroid_col"         : tables.UInt16Col(),
        
        "bounding_box_min_row" : tables.UInt16Col(),
        "bounding_box_min_col" : tables.UInt16Col(),
        "bounding_box_max_row" : tables.UInt16Col(),
        "bounding_box_max_col" : tables.UInt16Col()
    }
    
    # NOTE if using flat field correction, make these floats
    # Otherwise, integers
    for c in channels:
        # Total intensities
        column_types["total_intensity_{}".format(c)] = tables.UInt32Col()
    
    # Need to know the axmimum filename length to make a string column
    
    max_channel_name_length = get_max_length(seg_channels)
    column_types["info_seg_channel"] = tables.StringCol(max_channel_name_length)
    
    max_filename_length = get_max_length(file_names)
    column_types["info_filename"] = tables.StringCol(max_filename_length)
    
    return column_types

# Input a list, return the length of the longest element
def get_max_length(items):
    longest = 0
    
    for n in items:
        l = len(n) 
        if l > longest:
            longest = l
    
    return longest

# Input a list of pixel co-ordinates (as per skimage.measure.regionpprops),
# an associated image (2d numpy array) to read intensities from, and and a background signal value.
# Return the sum of the pixel intensities minus the product of the background value and the number of pixels.
# NOTE is it faster to use the coords[] array, or to "multiply" by the mask?
def subtract_background_from_coords(coords, image, background):
    summed_intensity = image[coords[0], coords[1]].sum()
    total_background = len(coords) * background
    
    return summed_intensity - total_background

# Input an image, a binary mask, and a background value per pixel
# Return the sum of the pixel intensities minus the product of the background value and the number of pixels.
def subtract_background_from_region(image, image_mask, background):
    # Use mask to zero out unwanted pixels, and sum the intensities of the remaining ones
    summed_intensity = (image * image_mask).sum()
    
    # Calculate background signal
    # Sum of a mask is equivalent to the number of valid pixels, if the mask is as integer and not bool.
    total_background = image_mask.astype(numpy.uint8).sum() * background
    
    return summed_intensity - total_background

# Merge the tables at the very end
def merge_tables(in_file, out_file, cell_dict):
    Cell = make_cell_type(cell_dict)
    
    h5file_in = tables.open_file(in_file, mode="r")
    
    h5file_out = tables.open_file(out_file, mode="w")
    big_table = h5file_out.create_table(h5file.root, "concatenated_measurements", Cell)
    
    ### See: https://stackoverflow.com/questions/24891014/append-all-rows-from-one-table-to-another-using-pytables
    for small_table in h5file_in.walk_nodes(where="/", classname="Table")
        for row in small_table.iterrows():
            big_table.append([row[:]]) # row.read() ?

    # Done copying trench properties
    big_table.flush()
    big_table.close()
    
    h5file_in.close()
    h5flie_out.close()
