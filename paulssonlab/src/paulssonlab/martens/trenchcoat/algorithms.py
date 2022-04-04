import numpy

from skimage.morphology import remove_small_objects, binary_closing, binary_erosion
from skimage.segmentation import watershed
from skimage.measure import label
from skimage.filters import gaussian, threshold_otsu, threshold_niblack

# from ipywidgets import fixed, interactive

from scipy import ndimage

# TODO: finish the basic dual or multi-channel thresholding algo

### Single threshold


def run_single_threshold(stack, ch_to_index, params):
    """
    Very simple max thresholding example.
    """
    data = stack[ch_to_index[params["channel"]]]
    return single_threshold(data, **params["parameters"])


# def run_single_threshold_interactive(stack, ch_to_index, params):
# """
# Very simple max thresholding example.
# Interactive version, with ipywidgets.
# """
# data = stack[ch_to_index[params["channel"]]]
# return interactive(single_threshold, data=fixed(data), **params["parameters"])


def single_threshold(data, cutoff):
    return data < cutoff


### Niblack on fluorescence, & thresholding on phase


def run_fluor_phase_segmentation(stack, ch_to_index, params):
    """
    A version which uses thresholding on fluor & phase channels
    Thresholding which also uses Niblack (mean, std dev) & watershedding
    """
    stack_fl = stack[..., ch_to_index[params["fluorescent_channel"]]]
    stack_ph = stack[..., ch_to_index[params["phase_channel"]]]

    return fluor_phase_segmentation(stack_fl, stack_ph, **params["parameters"])


# def run_fluor_phase_segmentation_interactive(stack, ch_to_index, params):
# """
# A version which uses thresholding on fluor & phase channels
# Thresholding which also uses Niblack (mean, std dev) & watershedding.
# Interactive version, with ipywidgets.
# """
# stack_fl = stack[..., ch_to_index[params["fluorescent_channel"]]]
# stack_ph = stack[..., ch_to_index[params["phase_channel"]]]

# return interactive(
# fluor_phase_segmentation,
# stack_fl=fixed(stack_fl),
# stack_ph=fixed(stack_ph),
# **params["parameters"],
# )


def fluor_phase_segmentation(
    stack_fl,
    stack_ph,
    otsu_multiplier,
    garbage_otsu_value,
    fluor_sigma,
    phase_sigma,
    phase_threshold_min,
    phase_threshold_max,
    scaling_factor,
    niblack_w,
    niblack_k,
    min_size,
    fluor_background,
    phase_background,
):
    """
    An algorithm which can use both fluorescence and transmitted light (e.g. phase contrast)
    information to detect & distinguish cells from a stack of images.

    The algorithm is written to support stacks of frames, which is useful for
    handling stacks of trenches from a same mother machine micrograph.

    Some of the steps, such as the blurring and the phase mask, can be ignored
    if the parameters are set accordingly.
    """
    # NOTE: worry about integer underflow?
    stack_fl -= fluor_background
    stack_ph -= phase_background

    # Pre-apply a slight Gaussian blur, to help when occasional pixels dip below the threshold
    blurred = numpy.empty(stack_fl.shape, dtype=numpy.float32, order="F")
    for z in range(blurred.shape[2]):
        blurred[..., z] = gaussian(
            stack_fl[..., z],
            sigma=fluor_sigma,
            mode="reflect",
            multichannel=False,
            preserve_range=True,
        )

    # Threshold is calculated using Otsu method, which samples the distribution
    # of pixel intensities, and then scaled by the otsu multiplier, which might
    # differ depending on empirical experience with a given channel or dataset.
    #
    # Then set pixels < threshold to zero, but preserve values > threshold.
    # Helps with the mean + stdev*k calculations in Niblack.
    # It is "bimodal" because the values are either zero, or they are the
    # original value.
    # NOTE: is the broadcasting operating as expected for threshold[z] vs bimodal[..., z]?
    # bimodal = ((stack_fl >= threshold) * stack_fl).astype(stack_fl.dtype, casting="safe", copy=False)
    bimodal = numpy.empty(stack_fl.shape, dtype=numpy.float32, order="F")
    for z in range(bimodal.shape[2]):
        threshold = threshold_otsu(blurred[..., z]) * otsu_multiplier

        bimodal[..., z] = ((stack_fl[..., z] >= threshold) * stack_fl[..., z]).astype(
            stack_fl.dtype, casting="safe", copy=False
        )

        # Without using an if statement, multiply: bimodal *= (threshold > garbage_otsu_value)
        # Either x 1, stays the same, or x 0, and then all subsequent steps do nothing.
        bimodal[..., z] *= threshold > garbage_otsu_value

    del blurred
    del threshold

    # A mask from the phase channel
    phase_blur = numpy.empty(stack_fl.shape, dtype=numpy.float32, order="F")
    for z in range(phase_blur.shape[2]):
        phase_blur[..., z] = gaussian(
            stack_ph[..., z],
            sigma=phase_sigma,
            mode="reflect",
            multichannel=False,
            preserve_range=True,
        )

    phase_mask = (phase_blur > phase_threshold_min) & (phase_blur < phase_threshold_max)

    del phase_blur

    # Apply the phase mask
    bimodal *= phase_mask

    del phase_mask

    # Apply Niblack method
    # The "multiply by scaling_factor" trick emphasizes the values of non-zero pixels.
    # This enhances the std. dev. at the edges of cells.
    niblack = numpy.empty(stack_fl.shape, dtype=numpy.float32, order="F")
    for z in range(niblack.shape[2]):
        niblack[..., z] = threshold_niblack(
            bimodal[..., z] * scaling_factor, window_size=niblack_w, k=niblack_k
        )

    # Only keep pixels which are < than the mean*stdev*k (in rectangular region, size w).
    # NOTE does it help to specify it explicitly, with F ordering?
    # or does it actually make things worse, because the next line then trashes it and makes a new array???
    bimodal = bimodal.astype(numpy.bool, casting="unsafe", copy=False)
    mask = (stack_fl > niblack) * bimodal

    # Fill in internal holes in the "mask"
    for z in range(mask.shape[2]):
        mask[..., z] = ndimage.binary_fill_holes(mask[..., z])

    del niblack
    del bimodal

    # The goal of these steps is to make a single, small region within each cell
    # and to use these regions as seeds for the watershed algorithm.
    markers = numpy.empty(stack_fl.shape, dtype=numpy.uint16, order="F")
    selem = numpy.array([[1, 0, 1], [1, 0, 1], [1, 0, 1]])
    for z in range(markers.shape[2]):
        # Use erosion to remove one perimeter's worth of pixels.
        # The goal here is to help separate adjacent cells which may
        # be joined by a single pixel.
        # TODO determine whether or not this is beneficial?
        # Does it sometimes cause cells to "disappear" if it fully
        # erodes the labeled seed pixels?
        m = binary_erosion(mask[..., z], selem=selem)

        # Remove small objects.
        # Connectivity of 1 seems to help a lot in preventing adjacent cells from being merged?
        # But connectivity of 2 can help if the eroded bits are very small.
        remove_small_objects(
            m,
            # mask[..., z],
            min_size=10,
            connectivity=2,
            in_place=True,
        )

        # Label the regions
        # NOTE is there a way to change the labeling so that it
        # coincides with the "top down" view that we associate with trenches?
        # (it seems to go more left-to-right, which doesn't make sense!)
        # Yes, by making an intermediate, rotated labels array, and then
        # transposing it back & storing it again. Could maybe improve
        # performance by tweaking their algo. to operate in the correct dimensional order?
        markers[..., z] = label(
            m.T,
            # mask[..., z].T,
            connectivity=2,
        ).T

    del selem

    # Make a basin (invert max, min pixels). Use the original intensity image for
    # the pixel values to use for the gradient ascent. Take advantage of the mask
    # to exclude undesired pixels.
    # NOTE it would be nice if the watershed algo. allowed for a max priority queue.
    # Would save us a step here.
    # NOTE it shouldn't matter if we use the max across all regions, or a same region,
    # it only matters that it's truly larger than all other values.
    image = stack_fl.max() - stack_fl

    # Store a stack of labeled segmentation masks
    # Must be integer type (what is returned by watershedding, i.e. integer labels.)
    result = numpy.empty(stack_fl.shape, dtype=numpy.uint16, order="F")
    for z in range(result.shape[2]):
        # FIXME this fixes really odd behavior where
        # the watershed refuses to handle the mask properly
        # unless I multiply it by a bunch of boolean ones.
        o = numpy.ones(mask[..., z].shape, dtype=numpy.bool)
        result[..., z] = watershed(
            image=image[..., z], markers=markers[..., z], mask=mask[..., z] * o
        )

    # TODO Remove small objects again?

    return result


### Niblack on fluorescence
def run_niblack_segmentation(stack, ch_to_index, params):
    """
    Thresholding, Niblack (mean, std dev) & watershedding on a single channel
    Input: 4-D numpy array (channel, stack, X, Y)
           dict converting channel names to channel position index in the stack
           segmentation params
    Output: stack of binary masks
    FIXME/TODO: some tweaks haven't made it over to this version yet, see version with phase info.
    """
    stack = stack[ch_to_index[params["channel"]]]

    # Store a stack of labeled segmentation masks
    result = numpy.empty(stack.shape, dtype=numpy.uint16)

    # For each image in the stack
    for z, stack_elem in enumerate(stack):
        # Threshold is calculated using Otsu method, which samples the distribution of pixel intensities, and then
        # scaled by the otsu multiplier, which might differ depending on empirical experience with a given channel or dataset.
        threshold = threshold_otsu(stack_elem) * params["otsu_multiplier"]

        # If it's a low value, then it's probably just noise...
        if threshold < params["garbage_otsu_value"]:
            result[z] = numpy.zeros(shape=(stack.shape[1], stack.shape[2]))
        else:
            # Set pixels < threshold to zero, but preserve values > threshold.
            # Helps with the mean + stdev*k calculations in Niblack.
            # It is "bimodal" because the values are either zero, or they are the original value.
            # The "multiply by scaling_factor" trick helps to emphasize the value of the pixels which are kept.
            # This enhances the std. dev. at the edges of cells.
            # FIXME what happens when multiplying by a floating point scaling factor?
            bimodal = (stack_elem < threshold) * 0 + (
                stack_elem >= threshold
            ) * stack_elem.astype(numpy.uint32) * params["scaling_factor"]

            # Apply Niblack method
            niblack = threshold_niblack(
                bimodal, window_size=params["niblack_w"], k=params["niblack_k"]
            )

            # Only keep pixels which are < than the mean*stdev*k (in rectangular region, size w).
            mask = stack_elem < niblack

            # Closing helps to fill in internal holes inside the cell
            mask = binary_closing(mask)

            # Remove small objects
            remove_small_objects(
                mask, min_size=params["min_size"], connectivity=2, in_place=True
            )

            # Label the regions
            markers = label(mask)

            # Make a basin (invert max, min pixels)
            # NOTE it would be nice if the watershed algo. allowed for a max priority queue.
            # Would save us a step here.
            image = stack_elem_fl.max() - stack_elem_fl

            result[z] = watershed(image=image, markers=markers, mask=mask)

    return result


def measure_whole_trench(stack, ch_to_index, params):
    """
    Returns a single labeled region for each whole trench (no cell segmentation).
    Useful for measuring total intensities within a trench.
    Doesn't require any parameters.
    """
    # Yes, this is a simple as returning all 1s!
    new_shape = (stack.shape[0], stack.shape[1], stack.shape[2])
    return numpy.ones(new_shape, dtype=numpy.uint16)
