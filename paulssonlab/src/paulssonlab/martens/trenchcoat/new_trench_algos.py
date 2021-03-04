#!/usr/bin/env python3
import tables
import numpy
from scipy.signal import find_peaks, savgol_filter
import os
import pathlib

# from skimage.filters import unsharp_mask

# FIXME something is wrong with padding ~ I think the negative padding has been double negatived??

# TODO:
# There is clearly an offset between the YFP & MCHERRY channels.
#
# When is the best time to apply the sorts of corrections?
#
# How to determine a correction in brightfield? Will this effect be visible?

# NOTE: it seems like curvature along the x dimension in an FOV
# is enough to warp the image such that the apparent spacing between
# adjacent trenches actually isn't constant.
# Is there a good way to generate a set of boxes,
# place it using the overall optimization, and then tweak each box by nudging it a bit?
# For example, could check +/- 1 or 2 pixels left or right and see if its edges can be further minimized.


def launch(file, frames, pool, metadata, params, in_file, out_dir, pbar):
    """
    This function is directly called by Trenchcoat's trench_detect routine.
    It's up to the trench algo for how to iterate etc., but we provide
    a process pool, detection params dict,
    """
    # Generate a set of rectangles which represent trenches
    # but which do not yet have known x,y starting points
    # which best match one or more of the trench rows.
    # Re-use the same boxes for all iterations!
    p = params["parameters"]
    boxes = generate_boxes(
        p["trenches"]["tr_width"],
        p["trenches"]["tr_height"],
        p["trenches"]["tr_spacing"],
        metadata["width"],
    )

    # For the progress bar
    def update_pbar(*a):
        pbar.update()

    for fov in metadata["fields_of_view"]:
        for frame in frames:
            for z in metadata["z_levels"]:
                args = [in_file, out_dir, file, fov, frame, z, params, boxes]
                pool.apply_async(run_trench_analysis, args, callback=update_pbar)
                # run_trench_analysis(*args) # DEBUG


def generate_boxes(tr_width, tr_height, tr_spacing, img_width):
    """
    Input trench width, trench height, trench spacing and image width.
    Return sets of [min_row, 0, max_row, tr_height],
    enough to fill the entire width of the image without running over the edge.
    Note that, for our purposes, the min_col & max_col are redundant.

    Can call this multiple times & aggregate the results into a higher-order list or array, if needed.
    """
    xs = []

    # Keep track of iterations to account for fraction pixel offsets
    iter_num = 0
    pos = 0
    while pos < img_width - tr_width - int(tr_spacing):
        pos = int((tr_spacing + tr_width) * iter_num)
        xs.append([pos, 0, pos + tr_width, tr_height])
        iter_num += 1

    return numpy.array(xs)


def run_trench_analysis(in_file, out_dir, filename, fov, frame, z_level, params, boxes):
    """
    Run basic trench analysis on an HDF5 file: detect trenches for each time
    frame, or just using the first time frame. Then, measure trench properties
    in all time frames.

    Store detected trenches in a PyTables array.

    Each array entry (row) has 4 values (min_row, min_col, max_row, max_col), and there is 1 row per region (trench).
    Regions must be rectangular, and must have the same dimension (for stacking purposes).
    """
    # 1. Load the image, optionally crop it
    img = load_image_for_regions_analysis(
        in_file, filename, fov, frame, z_level, params
    )

    # 2. Create an HDF5 file for writing trench row & trench regions info.
    dir_path = os.path.join(out_dir, "{}/FOV_{}/Frame_{}/".format(filename, fov, frame))
    pathlib.Path(dir_path).mkdir(parents=True, exist_ok=True)

    out_file_path = os.path.join(dir_path, "Z_{}_regions.h5".format(z_level))
    regions_file = tables.open_file(out_file_path, mode="w")

    # 3. Analyze the image to find trench rows & trenches, write coordinates of the detected rows & trenches to HDF5
    # This method detects both the rows & trenches in a single step
    # FIXME no correction for cropping; and the results are NOT written to disk!
    # To do so, need to retroactively parse out the row coords from the trenches.
    p = params["parameters"]
    (trenches, y_offsets) = find_trenches(img, boxes, **p["trenches"])

    # Take the y_offsets and create bounding boxes for the rows of trenches
    rows = []
    for y in y_offsets:
        bbox = [0, y, img.shape[0], y + p["trenches"]["tr_height"]]
        rows.append(bbox)

    # Write all row coords to disk as a single array
    regions_file.create_array("/", "row_coords", obj=rows)

    # Write out the trench coordinates using a variable-length array in HDF5
    results = regions_file.create_vlarray(
        "/",
        "trench_coords",
        atom=tables.UInt16Atom(shape=(4)),
        expectedrows=len(trenches),
    )

    # Iterate the rows
    for r in trenches:
        results.append(r)

    # Done!
    regions_file.close()


def load_image_for_regions_analysis(in_file, filename, fov, frame, z_level, params):
    """
    Open the image to be used for features detection, & optionally apply cropping.
    """
    h5file = tables.open_file(in_file, mode="r")
    node_path = "/Images/{}/FOV_{}/Frame_{}/Z_{}/{}".format(
        filename, fov, frame, z_level, params["channel"]
    )

    # NOTE image was written in F-order. When reading from disk, image appears transposed.
    # For best performance, it makes sense to not re-transpose the image,
    # but instead to accept this new reference frame & interpret left/right/top/bottom,
    # and axis 0,1 appropriately.

    if "crop" in params:
        try:
            img = h5file.get_node(node_path)[
                params["crop"]["top"] : params["crop"]["bottom"],
                params["crop"]["left"] : params["crop"]["right"],
            ]
        except:
            # TODO Error!
            print(
                "Error loading cropped image. Invalid node? Missing crop params? Check channel name?"
            )
            pass
    else:
        # Load the whole image
        # TODO is there a way to get it to read into an empty numpy array with F-ordering?
        # Otherwise, we either have to rotate the image, or flip all of the co-ordinates.
        img = h5file.get_node(node_path).read()

        # And define the missing crop params, so that they are always defined (for convenience)
        # FIXME width & height? 0 & 1, or 1 & 0?
        params["crop"] = {
            "top": 0,
            "bottom": img.shape[0],
            "left": 0,
            "right": img.shape[1],
        }

    h5file.close()

    return img


def find_trenches(
    img,
    boxes,
    tr_width,
    tr_height,
    tr_spacing,
    pad_left,
    pad_right,
    pad_top,
    pad_bottom,
):
    """
    Find the optimal place(s) in the image where trenches
    with the given dimensions (width, height, spacing) fit.

    Return 2 things:
    1. the padded trenches as a Python list of numpy arrays,
    2. the min y-coordinates for the region(s) containing said trenches.
    """
    # Exhaustively compute the "matches" between all valid y offsets,
    # and 2 trench width of possible x offsets.
    (sums, ranges_y, ranges_x) = sum_trench_edge_offsets(
        img, boxes, tr_width, tr_height, tr_spacing
    )

    # For every y range, calculate the greatest difference in column sums along y,
    # comparing each x.
    # Peaks in this "space" are expected to occur where trenches,
    # which alternate between bright & dark regions, are to be found.
    g_diffs = greatest_diffs(sums, ranges_y, ranges_x)
    y_offsets = get_row_peaks(g_diffs)

    # Find the x offset
    rows_of_trenches = find_best_x_offset(y_offsets, sums, boxes)

    # Apply padding & filter out trenches which are out of bounds
    padded = pad_and_remove_out_of_bounds_trenches(
        img, rows_of_trenches, pad_left, pad_right, pad_top, pad_bottom
    )

    return (padded, y_offsets)


def sum_trench_edge_offsets(img, boxes, tr_width, tr_height, tr_spacing, axis=1):
    """
    Input pre-computed "boxes" representing rectangular co-ordinates
    (min_row, min_col, max_row, max_col) of trenches.
    Return all possible (x  * y) combinations of sums of pixel intensities along
    the left & right edges.
    Also return the maximum y-value (ranges_y) and x-value (ranges_x).
    x represents all possible offsets within 1 range of widths.
    y represents all possible y values withouth going past the edge of the image.

    If we calculate the cumulative along an axis, then
    the value of [x] - [x-length-1] == the sum of the values from [x-length] to [x].
    This method permits us to dramatically reduce the number of calculations.
    """
    # Pre-compute many of the sums, since they will be repeated for every range of y values
    cumulative = numpy.cumsum(img, axis=axis)

    # Don't go over the edge of the image!
    ranges_y = img.shape[1] - tr_height

    # FIXME what's the best ranges_x?
    # Explore 1 trench width + 1 trench spacing; or maybe 2 trench widths?
    ranges_x = tr_width * 2

    # FIXME test truncating a bit in the y direction from both ends, to avoid
    # the bright spots at trench entrances in transmitted light images.
    truncate_y = 10

    # Filter out any trenches which are within the bounds of the image,
    # but which would extend beyond the image when performing the sums
    boxes_inbounds = boxes[boxes[..., 2] + ranges_x < img.shape[0]]

    sums = numpy.empty(shape=(ranges_x, ranges_y), dtype=numpy.uint32)
    # Iterate all possible trench row values (without running over the edge of the image)
    for y in range(ranges_y):
        # Iterate ranges_x worth of x values (possibly 2 trench widths?)
        for x in range(ranges_x):
            # Extract pairs of x ranges (the edges of the trenches),
            # and add to them the current x value
            tr = boxes_inbounds.T
            left = tr[0] + x
            right = tr[2] + x

            top = y - truncate_y + tr_height
            bottom = y + truncate_y

            # Take advantage of operations on slices: the sums will be applied across all
            # desired edge values stored in left and right.
            # Also take advantage of the pre-computed cumulative array to reduce
            # redundant calculations.
            sums[x, y] = (
                cumulative[left, top].sum()
                - cumulative[left, bottom].sum()
                + cumulative[right, top].sum()
                - cumulative[right, bottom].sum()
            )

    return (sums, ranges_y, ranges_x)


def greatest_diffs(sums, ranges_y, ranges_x):
    """
    Exhaustively search y ranges, and a subset of x ranges, to calculate
    the greatest difference between x offset values.
    The idea is that starting y values will define the upper edge of a trench
    row, and that trench rows, which alternate between bright and dark spots,
    will have very large greatest differences even when shifting just a few pixels in x.
    """
    greatest_differences = numpy.empty(shape=ranges_y, dtype=numpy.uint32)

    for y in range(ranges_y):
        sorted_sums = numpy.sort(sums[..., y])
        greatest_differences[y] = sorted_sums[-1] - sorted_sums[0]

    return greatest_differences


def get_row_peaks(g_diffs, prom_divisor=5, filter=False, window_length=11):
    """
    Use the "greatest differences" between the different possible x positions
    to determine where this difference is maximized.
    The notion is that areas outside of trenches will not have appreciable differences,
    whereas there will be large differences within a trench row.

    Prominence is used in the find_peaks algorithm to help weed out bogus peaks.
    By default, a value 10x less than the maximum value (greatest peak)
    is used to calculate this value.

    Optionally, the data can be run through a Savitzky-Golay smoothing filter.
    In principle this can help to smooth out local variations.
    In practice, it seems to introduce subtle biases which might
    improperly shift the peaks.
    """
    prom = g_diffs.max() / prom_divisor

    if filter:
        g_diffs = savgol_filter(g_diffs, window_length=window_length, polyorder=1)

    peaks, peaks_info = find_peaks(g_diffs, prominence=prom)

    return peaks


def find_best_x_offset(y_offsets, sums, boxes):
    """
    Now that we have valid peak(s), use those y-values to also pick
    the correct x-values.

    Each y offset corresponds to the y location of a peak ~ which should be
    the starting point for a trench row.

    Create a list. A set of trenches will be returned for every trench row.
    """
    result = numpy.empty(shape=(len(y_offsets), len(boxes), 4), dtype=numpy.uint16)
    # Iterate all trench rows
    for i, y_off in enumerate(y_offsets):
        # Which array index contains the smallest sum of
        # left, right border pixel values?
        # The index represents the amount to shift from zero.
        x_off = sums[..., y_off].argmin()

        # For each x & y offset, adjust all of the rectangles accordingly.
        result[i] = boxes + [x_off, y_off, x_off, y_off]

    return result


def rejiggle_trenches(img, globally_optimal):
    """
    Adjust each trench +/- 2 (TBD?) pixels (to the left or right)
    NOTE this doesn't actually seem to help much ~ might even make things worse.
    NOTE this function is currently not used.
    """
    # The range of x offsets to explore
    offsets = [i for i in range(-3, 4, 1)]

    truncate_y = 10

    # An array used to shift the left & right edges
    padding = numpy.zeros(shape=4, dtype=numpy.uint16)

    for row_num in range(globally_optimal.shape[0]):
        row = globally_optimal[row_num]
        for j, tr in enumerate(row):
            # All possible column sums for the offsets
            # NOTE is there a way to just map the trenches against the offsets?
            sums = numpy.empty(shape=len(offsets), dtype=numpy.uint16)
            for i, x in enumerate(offsets):
                s = (
                    img[tr[0] + x, tr[1] + truncate_y : tr[3] - truncate_y].sum()
                    + img[tr[2] + x, tr[1] + truncate_y : tr[3] - truncate_y].sum()
                )

                sums[i] = s

            # Find the x offset for the smallest sum, which
            # should align with the dark edges of a trench.
            x_off = sums.argmin()

            # Modify the padding for this trench
            padding = [offsets[x_off], 0, offsets[x_off], 0]

            # Apply the padding to the trench
            # FIXME what if this sends the trench out of bounds??
            row[j] += padding

    # Globally optimal is now rejiggled to be somewhat more locally optimal
    return globally_optimal


def pad_and_remove_out_of_bounds_trenches(
    img, rows_of_trenches, pad_left, pad_right, pad_top, pad_bottom
):
    """
    Input a 3D list: trench rows, trenches, (min_row, min_col, max_row_max_col),
    and padding amounts for left, right, top, bottom.
    Apply the padding & filter out trenches which are out of bounds.
    Return a list of numpy arrays.
    """
    # Note that overall result may be a ragged array, and will not neatly fit
    # into e.g. a numpy array, if any trenches were removed from some rows but
    # not others.
    result = []

    for row in rows_of_trenches:
        filtered_tr = []
        for tr in row:
            # Add padding
            a = [
                tr[0] - pad_left,
                tr[1] - pad_top,
                tr[2] + pad_right,
                tr[3] + pad_bottom,
            ]

            # Make sure none of the trenches is out of bounds
            if (
                (a[0] >= 0)
                & (a[0] < img.shape[0])
                & (a[1] >= 0)
                & (a[1] < img.shape[1])
                & (a[2] >= 0)
                & (a[2] < img.shape[0])
                & (a[3] >= 0)
                & (a[3] < img.shape[1])
            ):
                filtered_tr.append(a)

        result.append(numpy.array(filtered_tr, dtype=numpy.uint16))

    return result


# def find_trench_rows(img, boxes, tr_width, tr_height, tr_spacing):
# """
# These functions were originally intended to also find the trenches themselves.
# However, it is also possible to use them to find trench rows.
# Still, it's necessary to present some sort of trench geometry.
# TODO is it possible to further simplify this function?
# Could we just apply some sort of measure of variance within broad rectangular regions?
# NOTE this function is currently not used.
# """
## This seems to help enhance the edges around trenches
## TODO is it actually helpful enough to warrant including?
## What should the parameters be?
##img = unsharp_mask(img.astype(numpy.uint16), preserve_range=True)

## Exhaustively compute the "matches" between all valid y offsets,
## and 1 trench width of possible x offsets.
# (sums, ranges_y, ranges_x) = sum_trench_edge_offsets(img, boxes, tr_width, tr_height, tr_spacing)

## Calculate the x offsets which are most different for every y offset.
## Peaks in this "space" are expected to occur where trenches,
## which alternate between bright & dark regions, are to be found.
# g_diffs = greatest_diffs(sums, ranges_y, ranges_x)
# y_offsets = get_row_peaks(g_diffs)

# rows = numpy.empty(shape=(len(y_offsets), 4), dtype=numpy.uint16)

# for i, y in enumerate(y_offsets):
# rows[i] = [0, y, img.shape[0], y + tr_height]

# return rows
