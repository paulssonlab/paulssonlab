# Copying in old arrayfire code
# TODO: modify it to work with my new pipeline
# NOTE: there's an official padding suborutine now? double check
# https://github.com/arrayfire/arrayfire/pull/2682

# Set up the arrayfire environment
import arrayfire

# NOTE also need to pick a device
# arrayfire.set_backend('opencl') # TODO this has to be set in the main file?

# gauss_sigma = 0.5
# gauss_kernel = arrayfire.image.gaussian_kernel(3, 3, sigma_r=gauss_sigma, sigma_c=gauss_sigma)

## Convolution kernel
# ws = 11
# window_size = ws
# window_sq = ws ** 2
# kernel = make_niblack_kernel(ws, data_type)


### Segmentation operations


def load_stack(channel, stack, x_dimension, y_dimension, stack_size):
    """
    Load a stack of images of the same channel from a given image group, return an ArrayFire array
    """
    af_stack = arrayfire.array.Array(
        dims=(x_dimension, y_dimension, stack_size), dtype=arrayfire.Dtype.f32
    )
    for i, image_group in enumerate(stack):
        af_stack[:, :, i] = arrayfire.interop.from_ndarray(
            image_group._f_get_child(channel)[:]
        )

    return af_stack


def correct_stack(af_stack, channel, af_background_matrices):
    """
    Subtract camera noise and normalize by flat field correction
    """
    for z in arrayfire.ParallelRange(af_stack.shape[2]):
        af_stack[:, :, z] = (
            af_stack[:, :, z] - af_background_matrices["DARK"]
        ) / af_background_matrices[channel]

    return af_stack


def morphological_open(image):
    """
    Morphological open: erosion, then dilation
    """
    return arrayfire.image.dilate(arrayfire.image.erode(image))


def morphological_close(image):
    """
    Morphological close: dilation, then erosion
    """
    return arrayfire.image.erode(arrayfire.image.dilate(image))


def run_niblack_phase_segmentation_GPU(stack, ch_to_index, params):
    """
    Function to call segmentation.
    FIXME is this function defined twice????
    """
    stack_fl = stack[ch_to_index[params["fluorescent_channel"]]]
    stack_fl = arrayfire.interop.from_ndarray(stack_fl)

    stack_ph = stack[ch_to_index[params["phase_channel"]]]
    stack_ph = arrayfire.interop.from_ndarray(stack_ph)

    return niblack_phase_segmentation_GPU(stack_fl, stack_ph, **params["parameters"])


# C++ Code from http://arrayfire.org/docs/image_processing_2binary_thresholding_8cpp-example.htm:
# array otsu(const array& in) {
#    array gray;
#    int channels = in.dims(2);
#    if (channels > 1)
#        gray = colorSpace(in, AF_GRAY, AF_RGB);
#    else
#        gray = in;
#    unsigned total = gray.elements();
#    array hist     = histogram(gray, 256, 0.0f, 255.0f);
#    array wts      = range(256);
#
#    array wtB   = accum(hist);
#    array wtF   = total - wtB;
#    array sumB  = accum(wts * hist);
#    array meanB = sumB / wtB;
#    float lastElemInSumB;
#    sumB(seq(255, 255, 1)).host((void*)&lastElemInSumB);
#    array meanF = (lastElemInSumB - sumB) / wtF;
#    array mDiff = meanB - meanF;
#
#    array interClsVar = wtB * wtF * mDiff * mDiff;
#
#    float max        = af::max<float>(interClsVar);
#    float threshold2 = where(interClsVar == max).scalar<unsigned>();
#    array threshIdx  = where(interClsVar >= max);
#    float threshold1 =
#        threshIdx.elements() > 0 ? threshIdx.scalar<unsigned>() : 0.0f;
#
#    return threshold(gray, (threshold1 + threshold2) / 2.0f);
# }


def otsu_GPU(image):
    """
    Run Otsu thresholding on an image
    """
    total = image.elements()

    min_val = arrayfire.algorithm.min(image)
    max_val = arrayfire.algorithm.max(image)

    # FIXME: which one is right?
    nbins = int(max_val) - int(min_val)
    nbins = int(max_val - min_val)

    hist = arrayfire.image.histogram(image, nbins, min_val=min_val, max_val=max_val)
    weights = arrayfire.data.range(nbins)

    wt_b = arrayfire.algorithm.accum(hist)
    wt_f = total - wt_b

    sum_b = arrayfire.algorithm.accum(weights * hist)
    mean_b = sum_b / wt_b

    last_elem_in_sum_b = sum_b[-1]

    mean_f = (last_elem_in_sum_b - sum_b) / wt_f
    m_diff = mean_b - mean_f

    inter_cls_var = wt_b * wt_f * arrayfire.arith.pow(m_diff, 2)

    max = arrayfire.algorithm.max(inter_cls_var)
    threshold_2 = arrayfire.algorithm.where(inter_cls_var == max)[0]

    thresh_idx = arrayfire.algorithm.where(inter_cls_var >= max)
    if thresh_idx.elements() > 0:
        threshold_1 = thresh_idx[0]
    else:
        threshold_1 = 0.0

    return (threshold_1 + threshold_2) / 2.0


def get_phase_stack(
    stack,
    phase_channel_name,
    x_dimension,
    y_dimension,
    stack_size,
    af_background_matrices,
    gauss_kernel,
    phase_cutoff,
):
    """
    Load a stack of phase images (multiple image groups) & run a Gaussian blur + apply a cutoff.
    Return an ArrayFire array of phase images, from multiple image groups, which has been "binarized"
    depending on whether or not the pixels are greater than or less than the cutoff value.
    """
    stack_phase_af = load_stack(
        phase_channel_name, stack, x_dimension, y_dimension, stack_size
    )

    # FIXME best way to correct phase images for flat field?
    #     stack_phase_af = correct_stack(stack_phase_af, phase_channel_name, af_background_matrices)

    # Blur
    stack_phase_af = arrayfire.signal.convolve2(stack_phase_af, gauss_kernel)

    # Cutoff
    stack_phase_af = (stack_phase_af < phase_cutoff) * 1.0

    return stack_phase_af


# FIXME is this defined twice??
def run_niblack_phase_segmentation_GPU(stack, ch_to_index, params):
    stack_fl = stack[..., ch_to_index[params["fluorescent_channel"]]]
    stack_fl = arrayfire.interop.from_ndarray(stack_fl)
    stack_ph = stack[..., ch_to_index[params["phase_channel"]]]
    stack_ph = arrayfire.interop.from_ndarray(stack_ph)

    return niblack_phase_segmentation_GPU(stack_fl, stack_ph, **params["parameters"])


def niblack_phase_segmentation_GPU(
    stack_fl,
    stack_ph,
    garbage_otsu_value,
    phase_sigma,
    phase_threshold_min,
    phase_threshold_max,
    scaling_factor,
    niblack_w,
    niblack_k,
    fluor_background,
    phase_background,
    otsu_multiplier,
):
    """
    Run segmentation on the GPU using the ArrayFire.
    Input:
    1. a stack of corrected fluorescence images,
    2. a stack of corrected & pre-processed phase images,
    3. and some segmentation parameters
    Return 3 Numpy arrays necessary for the skimage watershed algorithm.
    TODO pass in the convolution kernel, or re-initialize it every time?
    """
    # Can we make the kernels global & avoid having to pass them in somehow? Avoid re-calculating them every time?

    # Simple background correction
    # NOTE: worry about integer underflow?
    stack_fl -= fluor_background
    stack_ph -= phase_background

    # Calculate a threshold using the statistics of the image
    # NOTE correct dtype?
    num_regions = stack_fl.shape[-1]
    threshold = arrayfire.array.Array(dims=(num_regions), dtype=arrayfire.Dtype.u16)
    for z in arrayfire.ParallelRange(num_regions):
        threshold[z] = otsu_GPU(stack_fl[..., z])

    threshold *= otsu_multiplier

    # TODO: test against garbage otsu value? return zeros if fails the test!

    # Set pixels < threshold to zero, but preserve values > threshold.
    # Helps with the mean + stdev*k calculations in Niblack.
    # It is "bimodal" because the values are either zero, or they are the original value.
    # NOTE how is broadcasting going to work here? 1 threshold per region...
    bimodal = ((stack_fl >= threshold) * stack_fl).as_type(stack_fl.dtype)

    # A mask from the phase channel
    gauss_kernel = arrayfire.image.gaussian_kernel(
        3, 3, sigma_r=gauss_sigma, sigma_c=gauss_sigma
    )
    phase_blur = arrayfire.signal.convolve2(stack_ph, gauss_kernel)
    phase_mask = (phase_blur > phase_threshold_min) & (phase_blur < phase_threshold_max)

    # Apply the phase mask
    bimodal *= phase_mask

    # Binary mask
    binary_bimodal = bimodal_stack > 0.0

    # Run an opening (erosion followed by dilation) to get rid of lone pixels and to smooth the cell contours
    binary_bimodal = morphological_open(binary_bimodal)

    # And apply this operator's result back onto the bimodal_stack
    bimodal_stack = (binary_bimodal == 0) * 0.0 + (binary_bimodal == 1) * bimodal_stack

    # Niblack method: only keep those pixels which are greater in the original image
    # and were "cellular" in the (bimodal fluorescence * phase) mask
    niblack_kernel = make_niblack_kernel(niblack_w, data_type)
    mask = (
        stack_fl
        > niblack_GPU(bimodal_stack, kernel, niblack_w, niblack_k, stack_fl.dtype)
    ) * binary_bimodal

    del bimodal_stack

    # Find "connected components" (i.e. cells) regions within each image region (i.e. trench)
    # by labeling the binary mask with unique identifying labels.
    # Must explicitly iterate, because we are not interesting in multi-dimensional volumes!
    # FIXME for some reason the parallel iterator does not work!
    # Is the slicing/indexing function broken?
    # for z in arrayfire.ParallelRange(num_regions):
    markers = arrayfire.array.Array(dims=mask.dims(), dtype=mask.dtype)
    for z in range(num_regions):
        markers[..., z] = arrayfire.image.regions(mask[..., z])

    mask = binary_bimodal.to_ndarray()
    markers = markers.to_ndarray()

    del binary_bimodal
    arrayfire.device.device_gc()

    # FIXME change this to subtracting the maximum from all pixels!
    # Multiply by -1 to make it a basin
    # And cast to signed integer because the skimage watershed algorithm requires it
    # basin = (stack_fl * -1).as_type(arrayfire.Dtype.s16).to_ndarray()
    image = arrayfire.algorithms.max(stack_fl) - stack_fl
    image = image.to_ndarray()

    # Store a stack of labeled segmentation masks
    # Must be integer type (what is returned by watershedding, i.e. integer labels.)
    result = numpy.empty(image.shape, dtype=numpy.uint16, order="F")
    for z in range(image.shape[2]):
        result[..., z] = skimage.segmentation.watershed(
            image=image[..., z], markers=markers[..., z], mask=mask[..., z]
        )

    return result


def get_regions(
    image_stack, binary_bimodal, bimodal_stack, kernel, window_size, k, dtype
):
    """ """
    # Assume that there are at least 2 images in a stack, therefore dims() returns a tuple size 3
    num_images = image_stack.dims()[2]

    # Niblack method
    # Only keep those pixels which are greater in the original image
    # and were "cellular" in the (bimodal fluorescence * phase) mask
    result = compare_niblack(
        image_stack, binary_bimodal, bimodal_stack, kernel, window_size, k, dtype
    )

    # Must explicitly iterate, because we are not interesting in multi-dimensional volumes!
    #     for z in arrayfire.ParallelRange(num_images):

    # FIXME for some reason the parallel iterator does not work!
    # Is the slicing/indexing function broken?
    for z in range(num_images):
        result[..., z] = arrayfire.image.regions(result[..., z])

    result = result.to_ndarray()

    # Freedom!
    arrayfire.device.device_gc()

    return result


def niblack_GPU(image_stack, kernel, window_size, k, dtype):
    """
    Niblack: mean - k * std. dev
    """
    mean = mean_GPU(image_stack, kernel, window_size, dtype)

    # Std. dev
    # First, calculate the difference of the sum of the squares from the square of the sum
    # Don't store the mean squared in a separate variable
    s = mean_GPU(
        arrayfire.arith.pow(image_stack, 2), kernel, window_size, dtype
    ) - arrayfire.arith.pow(mean, 2)

    # Freedom!
    arrayfire.device.device_gc()

    # Negative numbers may arise due to floating point error -> set to 0.0
    # Then, square root
    s = arrayfire.arith.sqrt((s <= 0.0) * 0.0 + (s > 0.0) * s)

    # Apply Niblack method
    result = mean - k * s

    # Freedom!
    del mean
    del s
    arrayfire.device.device_gc()

    return result


def mean_GPU(image_stack, kernel, window_size, dtype):
    """
    Mean pixel value within a square region
    """
    (y_dimension, x_dimension, num_images) = image_stack.dims()

    # Padding
    result = arrayfire.array.Array(
        dims=(y_dimension + window_size * 2, x_dimension + window_size * 2, num_images),
        dtype=dtype,
    )
    for z in arrayfire.ParallelRange(num_images):
        result = pad_af(result[..., z], image_stack[..., z], window_size, z)

    # Integral image
    result = arrayfire.image.sat(result)

    # Convolution
    result = arrayfire.signal.convolve2(result, kernel)

    # Cropping and dividing by the size of the window (square region)
    result = (
        result[
            window_size - 1 : -1 * window_size - 1,
            window_size - 1 : -1 * window_size - 1,
        ]
        / window_sq
    )

    return result


def pad_af(arr_new, arr_old, pad, dim):
    """
    Pad an array (useful for handling image edges during matrix convolutions)

    Allocate a single larger array, and then perform 8 copies:
    1 - 4: the top, bottom, left, and right -- each one of these flipped once
    and 5 - 8. the 4 corners, made by flipping one of the sides again.
    T = top; B = bottom; M = middle; L = left; R = right
    TL  TM  TR
    ML  MM  MR
    BL  BM  BR

    NOTE: this will not work if the window size is bigger than the array.
    Since we typically use w = 11, and array = 2048x2048, this doesn't matter.

    1. Array 2. padding width (the same for all 4 edges)
    """
    # FIXME even if passing a 2D slice, it returns 3 dimensions
    # Need to check whether 2 or 3
    (in_dim_y, in_dim_x, in_dim_z) = arr_old.dims()

    # NOTE: array[y1 : y2, x1 : x2] := y is top -> bottom, x is left -> right

    # For counting backwards from an edge
    margin = -1 * pad

    # Copy the inside
    arr_new[pad : pad + in_dim_y, pad : pad + in_dim_x, dim] = arr_old[:, :, dim]

    # MR
    arr_new[pad:margin, margin:] = arrayfire.flip(
        arr_old[:, margin - 1 : -1, dim], dim=1
    )

    # ML
    arr_new[pad:margin, :pad] = arrayfire.flip(arr_old[:, 1 : pad + 1, dim], dim=1)

    # BM
    arr_new[margin:, pad:margin] = arrayfire.flip(arr_old[:pad, :, dim], dim=0)

    # TM
    arr_new[:pad, pad:margin] = arrayfire.flip(arr_old[margin:, :, dim], dim=0)

    # TL
    arr_new[:pad, :pad] = arrayfire.flip(arr_new[:pad, pad : pad * 2, dim], dim=1)

    # TR
    arr_new[:pad, margin:] = arrayfire.flip(
        arr_new[:pad, margin * 2 : margin, dim], dim=1
    )

    # BL
    arr_new[margin:, :pad] = arrayfire.flip(arr_new[margin:, pad : pad * 2, dim], dim=1)

    # BR
    arr_new[margin:, margin:] = arrayfire.flip(
        arr_new[margin:, margin * 2 : margin, dim], dim=1
    )

    return arr_new


def make_niblack_kernel(window_size, dtype):
    """
    Make a convolution kernel for the ArrayFire Niblack algorithm using the given window size
    """
    kernel = arrayfire.data.constant(
        val=0.0, d0=window_size + 1, d1=window_size + 1, dtype=dtype
    )
    kernel[0, 0] = 1.0
    kernel[-1, 0] = -1.0
    kernel[0, -1] = -1.0
    kernel[-1, -1] = 1.0

    return kernel
