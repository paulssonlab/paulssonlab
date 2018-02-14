import numpy as np
import skimage


def hessian_eigenvalues(img):
    I = skimage.filters.gaussian(img, 1.5)
    I_x = skimage.filters.sobel_h(I)
    I_y = skimage.filters.sobel_v(I)
    I_xx = skimage.filters.sobel_h(I_x)
    I_xy = skimage.filters.sobel_v(I_x)
    I_yx = skimage.filters.sobel_h(I_y)
    I_yy = skimage.filters.sobel_v(I_y)
    kappa_1 = (I_xx + I_yy) / 2
    kappa_2 = (np.sqrt((I_xx + I_yy) ** 2 - 4 * (I_xx * I_yy - I_xy * I_yx))) / 2
    k1 = kappa_1 + kappa_2
    k2 = kappa_1 - kappa_2
    k1[np.isnan(k1)] = 0
    k2[np.isnan(k2)] = 0
    return k1, k2


def extract_kymograph(img_series, x0, x1):
    num_timepoints = img_series.shape[0]
    xs, ys = coords_along(x0, x1)
    kymo = np.zeros((len(xs), num_timepoints))
    for t in range(num_timepoints):
        kymo[:, t] = img_series[t, ys, xs][::-1]
    return kymo


def get_image_series_segmentation(img_series, x0, x1):
    pass  # list of trenches, for each trench, a list of cell masks


def map_over_segmentation(img_series, cell_seg, func):
    pass


def bounding_box(points):
    upper_left_x = min(point[0] for point in points)
    upper_left_y = min(point[1] for point in points)
    lower_right_x = max(point[0] for point in points)
    lower_right_y = max(point[1] for point in points)
    return (
        np.array([upper_left_x, upper_left_y]),
        np.array([lower_right_x, lower_right_y]),
    )


def crop_point(x, x_lim, y_lim):
    return np.clip(x, *np.vstack([[x_lim], [y_lim]]).T)


def get_trench_thumbnail(img, trench_points, trench_idx):
    x_lim, y_lim = get_img_limits(img)
    ul, lr = get_trench_bbox(trench_points, trench_idx, x_lim, y_lim)
    return img[ul[1] : lr[1], ul[0] : lr[0]]


def get_trench_bbox(trench_points, trench_idx, x_lim, y_lim):
    # trench_points[0][trench_idx], trench_points[1][trench_idx]
    num_trenches = min(len(trench_points[0]), len(trench_points[1]))
    if not 0 <= trench_idx < num_trenches:
        raise ValueError("trench index out of bounds")
    points = [trench_points[i][trench_idx] for i in (0, 1)]
    if trench_idx == 0:
        x0_prev = 2 * trench_points[0][0] - trench_points[0][1]
        x1_prev = 2 * trench_points[1][0] - trench_points[1][1]
        x0_prev = crop_point(x0_prev, x_lim, y_lim)
        x1_prev = crop_point(x1_prev, x_lim, y_lim)
        points += [x0_prev, x1_prev]
    else:
        points += [trench_points[i][trench_idx - 1] for i in (0, 1)]
    if trench_idx + 1 == num_trenches:
        x0_next = 2 * trench_points[0][-1] - trench_points[0][-2]
        x1_next = 2 * trench_points[1][-1] - trench_points[1][-2]
        x0_next = crop_point(x0_next, x_lim, y_lim)
        x1_next = crop_point(x1_next, x_lim, y_lim)
        points += [x0_next, x1_next]
    else:
        points += [trench_points[i][trench_idx + 1] for i in (0, 1)]
    return bounding_box(points)


def image_sharpness(img):
    # FROM: https://stackoverflow.com/questions/7765810/is-there-a-way-to-detect-if-an-image-is-blurry/7767755#7767755
    img_blurred = skimage.filters.gaussian(img, 1)
    img_lofg = skimage.filters.laplace(img_blurred)
    return np.percentile(img_lofg, 99.9)
