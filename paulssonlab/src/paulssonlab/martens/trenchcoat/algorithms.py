import numpy

import skimage.morphology
import skimage.segmentation
import skimage.measure
import skimage.filters

from ipywidgets import fixed, interactive

# TODO: finish the basic dual or multi-channel thresholding algo


def run_single_threshold(stack, ch_to_index, params):
    """
    Very simple max thresholding example.
    """
    data = stack[ch_to_index[params["channel"]]]
    return single_threshold(data, **params["parameters"])


def run_single_threshold_interactive(stack, ch_to_index, params):
    """
    Very simple max thresholding example.
    Interactive version, with ipywidgets.
    """
    data = stack[ch_to_index[params["channel"]]]
    return interactive(single_threshold, data=fixed(data), **params["parameters"])


def single_threshold(data, cutoff):
    return data < cutoff


def run_niblack_phase_segmentation(stack, ch_to_index, params):
    """
    A version which uses thresholding on fluor & phase channels
    Thresholding which also uses Niblack (mean, std dev) & watershedding
    """
    stack_fl = stack[..., ch_to_index[params["fluorescent_channel"]]]
    stack_ph = stack[..., ch_to_index[params["phase_channel"]]]

    return niblack_phase_segmentation(stack_fl, stack_ph, **params["parameters"])


def run_niblack_phase_segmentation_interactive(stack, ch_to_index, params):
    """
    A version which uses thresholding on fluor & phase channels
    Thresholding which also uses Niblack (mean, std dev) & watershedding.
    Interactive version, with ipywidgets.
    """
    stack_fl = stack[..., ch_to_index[params["fluorescent_channel"]]]
    stack_ph = stack[..., ch_to_index[params["phase_channel"]]]

    return interactive(
        niblack_phase_segmentation,
        stack_fl=fixed(stack_fl),
        stack_ph=fixed(stack_ph),
        **params["parameters"],
    )


def niblack_phase_segmentation(
    stack_fl,
    stack_ph,
    otsu_multiplier,
    garbage_otsu_value,
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
    TODO: I bet this would be much faster if I loop individual steps which much be looped,
    and use UFUNCS for all the other steps!

    NOTE: preliminary benchmarking suggests this version is actually a bit slower?
    BUT there is only 1 region, so this maybe isn't the best test case...

    AND maybe it'll expand better to ArrayFire...
    Also, didn't check the outputs to make sure they're identical.
    TODO: free up memory as soon as possible, to reduce overhead?
    """
    # NOTE: worry about integer underflow?
    stack_fl -= fluor_background
    stack_ph -= phase_background

    # Threshold is calculated using Otsu method, which samples the distribution of pixel intensities, and then
    # scaled by the otsu multiplier, which might differ depending on empirical experience with a given channel or dataset.
    threshold = numpy.empty(stack_fl.shape[2], dtype=numpy.float32, order="F")
    for z in range(stack_fl.shape[2]):
        threshold[z] = (
            skimage.filters.threshold_otsu(stack_fl[..., z]) * otsu_multiplier
        )

    # Set pixels < threshold to zero, but preserve values > threshold.
    # Helps with the mean + stdev*k calculations in Niblack.
    # It is "bimodal" because the values are either zero, or they are the original value.
    # NOTE: is the broadcasting operating as expected for threshold[z] vs bimodal[..., z]?
    # bimodal = ((stack_fl >= threshold) * stack_fl).astype(stack_fl.dtype, casting="safe", copy=False)

    bimodal = numpy.empty(stack_fl.shape, dtype=numpy.float32, order="F")
    for z in range(stack_fl.shape[2]):
        bimodal[..., z] = (
            (stack_fl[..., z] >= threshold[z]) * stack_fl[..., z]
        ).astype(stack_fl.dtype, casting="safe", copy=False)
        # TODO: verify that this works
        # Without using an if statement, is the best way to ... ? multiply: bimodal *= (threshold > garbage_otsu_value)
        # Either x 1, stays the same, or x 0, and then all subsequent steps do nothing.
        bimodal[..., z] *= threshold[z] > garbage_otsu_value

    del threshold

    # A mask from the phase channel
    phase_blur = numpy.empty(stack_fl.shape, dtype=numpy.float32, order="F")
    for z in range(stack_fl.shape[2]):
        phase_blur[..., z] = skimage.filters.gaussian(
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
    for z in range(stack_fl.shape[2]):
        niblack[..., z] = skimage.filters.threshold_niblack(
            bimodal[..., z] * scaling_factor, window_size=niblack_w, k=niblack_k
        )

    # Only keep pixels which are < than the mean*stdev*k (in rectangular region, size w).
    # NOTE does it help to specify it explicitly, with F ordering?
    # or does it actually make things worse, because the next line then trashes it and makes a new array???
    bimodal = bimodal.astype(numpy.bool, casting="unsafe", copy=False)
    mask = (stack_fl > niblack) * bimodal

    del niblack
    del bimodal

    markers = numpy.empty(stack_fl.shape, dtype=numpy.float32, order="F")
    for z in range(stack_fl.shape[2]):
        # Closing fills small holes
        mask[..., z] = skimage.morphology.binary_closing(mask[..., z])

        # Remove small objects
        skimage.morphology.remove_small_objects(
            mask[..., z], min_size=min_size, connectivity=2, in_place=True
        )

        # Label the regions
        markers[..., z] = skimage.measure.label(mask[..., z])

    # Make a basin (invert max, min pixels)
    # NOTE it would be nice if the watershed algo. allowed for a max priority queue.
    # Would save us a step here.
    # NOTE it shouldn't matter if we use the max across all regions, or a same region, it only matters that it's
    # truly larger than all other values.
    image = stack_fl.max() - stack_fl

    # Store a stack of labeled segmentation masks
    # Must be integer type (what is returned by watershedding, i.e. integer labels.)
    result = numpy.empty(stack_fl.shape, dtype=numpy.uint16, order="F")
    for z in range(stack_fl.shape[2]):
        result[..., z] = skimage.segmentation.watershed(
            image=image[..., z], markers=markers[..., z], mask=mask[..., z]
        )

    return result


def niblack_phase_segmentation_old(
    stack_fl,
    stack_ph,
    otsu_multiplier,
    garbage_otsu_value,
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
    NOTE this version is ever so slightly faster... almost seems by a constant ~0.5-1.0 secs,
    regardless of the # of processes?
    """

    ## DEBUG, just to see how quickly the code runs when doing a minimal amount of work
    # return numpy.zeros(stack_fl.shape, dtype=numpy.uint16, order='F')

    # Store a stack of labeled segmentation masks
    # Must be integer type (what is returned by watershedding, i.e. integer labels.)
    result = numpy.empty(stack_fl.shape, dtype=numpy.uint16, order="F")

    # NOTE: worry about integer underflow?
    stack_fl -= fluor_background
    stack_ph -= phase_background

    # For each image in the stack, (z usually represents a *region*, but we call it z for generic purposes)
    for z in range(stack_fl.shape[2]):
        stack_elem_fl = stack_fl[..., z]
        stack_elem_ph = stack_ph[..., z]

        # Threshold is calculated using Otsu method, which samples the distribution of pixel intensities, and then
        # scaled by the otsu multiplier, which might differ depending on empirical experience with a given channel or dataset.
        threshold = skimage.filters.threshold_otsu(stack_elem_fl) * otsu_multiplier

        # If it's a low value, then it's probably just noise...
        if threshold < garbage_otsu_value:
            result[..., z] = numpy.zeros(
                shape=stack_elem_fl.shape, dtype=numpy.uint16, order="F"
            )
        else:
            # Set pixels < threshold to zero, but preserve values > threshold.
            # Helps with the mean + stdev*k calculations in Niblack.
            # It is "bimodal" because the values are either zero, or they are the original value.
            bimodal = ((stack_elem_fl >= threshold) * stack_elem_fl).astype(
                stack_fl.dtype, casting="safe", copy=False
            )
            # NOTE Does copy flag actually make a difference?

            # A mask from the phase channel
            phase_blur = skimage.filters.gaussian(
                stack_elem_ph,
                sigma=phase_sigma,
                mode="reflect",
                multichannel=False,
                preserve_range=True,
            )
            phase_mask = (phase_blur > phase_threshold_min) & (
                phase_blur < phase_threshold_max
            )

            # Apply the phase mask
            bimodal *= phase_mask

            # Apply Niblack method
            # The "multiply by scaling_factor" trick emphasizes the values of non-zero pixels.
            # This enhances the std. dev. at the edges of cells.
            niblack = skimage.filters.threshold_niblack(
                bimodal * scaling_factor, window_size=niblack_w, k=niblack_k
            )

            # Only keep pixels which are < than the mean*stdev*k (in rectangular region, size w).
            mask = (stack_elem_fl > niblack) * bimodal.astype(
                numpy.bool, casting="unsafe", copy=True
            )

            # Closing fills small holes
            mask = skimage.morphology.binary_closing(mask)

            # Remove small objects
            skimage.morphology.remove_small_objects(
                mask, min_size=min_size, connectivity=2, in_place=True
            )

            # Label the regions
            markers = skimage.measure.label(mask)

            # Make a basin (invert max, min pixels)
            # NOTE it would be nice if the watershed algo. allowed for a max priority queue.
            # Would save us a step here.
            # BROADCASTABLE ?
            image = stack_elem_fl.max() - stack_elem_fl
            # FIXME for F order
            result[..., z] = skimage.segmentation.watershed(
                image=image, markers=markers, mask=mask
            )

    return result


###


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
        threshold = (
            skimage.filters.threshold_otsu(stack_elem) * params["otsu_multiplier"]
        )

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
            niblack = skimage.filters.threshold_niblack(
                bimodal, window_size=params["niblack_w"], k=params["niblack_k"]
            )

            # Only keep pixels which are < than the mean*stdev*k (in rectangular region, size w).
            mask = stack_elem < niblack

            # Closing helps to fill in internal holes inside the cell
            mask = skimage.morphology.binary_closing(mask)

            # Remove small objects
            skimage.morphology.remove_small_objects(
                mask, min_size=params["min_size"], connectivity=2, in_place=True
            )

            # Label the regions
            markers = skimage.measure.label(mask)

            # Make a basin (invert max, min pixels)
            # NOTE it would be nice if the watershed algo. allowed for a max priority queue.
            # Would save us a step here.
            image = stack_elem_fl.max() - stack_elem_fl

            result[z] = skimage.segmentation.watershed(
                image=image, markers=markers, mask=mask
            )

    return result
