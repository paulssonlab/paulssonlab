import numpy
from scipy.ndimage import binary_fill_holes
from skimage.filters import threshold_niblack, threshold_otsu, unsharp_mask
from skimage.measure import label
from skimage.morphology import remove_small_objects
from skimage.segmentation import watershed


def run_fluor_sharp_segmentation(stack, ch_to_index, params):
    """
    A version which uses thresholding on fluor & phase channels
    Thresholding which also uses Niblack (mean, std dev) & watershedding
    """
    stack_fl = stack[..., ch_to_index[params["fluorescent_channel"]]]

    return fluor_sharpen_segmentation(stack_fl, **params["parameters"])


def fluor_sharpen_segmentation(
    stack_fl,
    otsu_multiplier,
    garbage_otsu_value,
    niblack_w,
    niblack_k,
    min_size,
    unsharp_mask_radius,
    unsharp_mask_amount,
):
    """
    Use fluorescence.
    Key difference from other Niblack-based routines is the use of sharpening.
    """
    # Sharpen the fluorescence image
    for z in range(stack_fl.shape[2]):
        stack_fl[..., z] = unsharp_mask(
            stack_fl[..., z],
            radius=unsharp_mask_radius,
            amount=unsharp_mask_amount,
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
        threshold = threshold_otsu(stack_fl[..., z]) * otsu_multiplier

        bimodal[..., z] = ((stack_fl[..., z] >= threshold) * stack_fl[..., z]).astype(
            stack_fl.dtype, casting="safe", copy=False
        )

        # Without using an if statement, multiply: bimodal *= (threshold > garbage_otsu_value)
        # Either x 1, stays the same, or x 0, and then all subsequent steps do nothing.
        bimodal[..., z] *= threshold > garbage_otsu_value

    # Apply Niblack method
    niblack = numpy.empty(stack_fl.shape, dtype=numpy.float32, order="F")
    for z in range(niblack.shape[2]):
        niblack[..., z] = threshold_niblack(
            bimodal[..., z], window_size=niblack_w, k=niblack_k
        )

    # Only keep pixels which are < than the mean*stdev*k (in rectangular region, size w).
    # NOTE does it help to specify it explicitly, with F ordering?
    # or does it actually make things worse, because the next line then trashes it and makes a new array???
    bimodal = bimodal.astype(numpy.bool, casting="unsafe", copy=False)
    mask = (stack_fl > niblack) * bimodal

    markers = numpy.empty(stack_fl.shape, dtype=numpy.uint16, order="F")
    for z in range(mask.shape[2]):
        # Fill in internal holes in the "mask"
        mask[..., z] = binary_fill_holes(mask[..., z])

        # Remove small objects.
        # Connectivity of 1 seems to help a lot in preventing adjacent cells from being merged?
        # But connectivity of 2 can help if the eroded bits are very small.
        mask[..., z] = remove_small_objects(
            mask[..., z], min_size=min_size, connectivity=1
        )

        # Label the regions
        # NOTE is there a way to change the labeling so that it
        # coincides with the "top down" view that we associate with trenches?
        # (it seems to go more left-to-right, which doesn't make sense!)
        # Yes, by making an intermediate, tranposed labels array, and then
        # transposing it back & storing it again. Could maybe improve
        # performance by tweaking their algo. to operate in the correct dimensional order?
        markers[..., z] = label(mask[..., z].T, connectivity=1).T

    del niblack
    del bimodal

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
            image=image[..., z],
            markers=markers[..., z],
            mask=mask[..., z] * o,  # multiply by boolean 1s
        )

    return result
