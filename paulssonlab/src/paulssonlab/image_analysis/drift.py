import warnings


def trench_cell_endpoints(img, sigma=2, k=2, min_height=0.3, margin_factor=1):
    img = skimage.img_as_float(img)
    profile = img.mean(axis=1)
    grad = misc.holoborodko_diff.holo_diff(
        1, scipy.ndimage.gaussian_filter1d(profile, sigma)
    )
    with warnings.catch_warnings(
        action="ignore", category=scipy.signal._peak_finding_utils.PeakPropertyWarning
    ):
        pos_peaks, pos_peak_props = scipy.signal.find_peaks(
            grad, height=min_height * grad.max(), width=(None, None)
        )
        neg_peaks, neg_peak_props = scipy.signal.find_peaks(
            -grad, height=-min_height * grad.min(), width=(None, None)
        )
    if len(pos_peaks) < 2:
        return None
    y1 = pos_peaks[0]
    y2 = neg_peaks[-1]
    margin = int(
        np.ceil(
            margin_factor
            * (pos_peak_props["widths"][0] + neg_peak_props["widths"][-1])
            / 2
        )
    )
    cutoff1 = min(y1 + 1 + margin, img.shape[0])
    cutoff2 = max(y2 - margin, 0)
    x1 = img[:cutoff1, :].mean(axis=0).argmax()
    x2 = img[cutoff2:, :].mean(axis=0).argmax()
    return np.array([[x1, y1], [x2, y2]])


class TranslationTransform(skimage.transform.EuclideanTransform):
    def estimate(self, src, dst):
        translation = (dst - src).mean(axis=0)
        self.params[0 : self.dimensionality, self.dimensionality] = translation
        return True


def ransac_translation(
    data, residual_threshold=3, min_samples=10, diagnostics=None, **kwargs
):
    data = (*np.array(data).swapaxes(0, 1),)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "divide by zero encountered in scalar divide")
        model_robust, inliers = skimage.measure.ransac(
            data,
            TranslationTransform,
            min_samples=min_samples,
            residual_threshold=residual_threshold,
            **kwargs,
        )
    if diagnostics is not None:
        feature_points = [
            [(*d1, "red" if inlier else "gray"), (*d2, "red" if inlier else "gray")]
            for d1, d2, inlier in zip(data[0], data[1], inliers)
        ]
        diagnostics["correspondences"] = hv.Path(feature_points, vdims=["color"]).opts(
            color="color", line_width=2
        )
    return model_robust.translation


def find_feature_drift(
    img1,
    img2,
    trenches,
    initial_shift=None,
    tolerance=1,
    feature_func=trench_cell_endpoints,
    estimation_func=ransac_translation,
    max_iterations=1,
    diagnostics=None,
):
    image_limits = geometry.get_image_limits(img1.shape)
    features1 = {}
    if initial_shift is None:
        initial_shift = np.array([0, 0])
    shifted_trenches = filter_trenches(
        shift_trenches(trenches, initial_shift), image_limits
    )
    for roi_idx, crop, ul in geometry.iter_crops(img1, shifted_trenches, corner=True):
        if (feature := feature_func(crop)) is not None:
            features1[roi_idx] = feature + ul[np.newaxis, ...]
    shift = initial_shift
    for i in range(max_iterations):
        features_list = []
        features2 = {}
        shifted_trenches = filter_trenches(
            shift_trenches(trenches, shift), image_limits
        )
        for roi_idx, crop, ul in geometry.iter_crops(
            img2, shifted_trenches, corner=True
        ):
            if (feature := feature_func(crop)) is not None:
                features2[roi_idx] = feature + ul[np.newaxis, ...]
        for roi_idx in features1.keys() & features2.keys():
            roi_features1 = features1[roi_idx]
            roi_features2 = features2[roi_idx]
            if roi_features1 is None or roi_features2 is None:
                continue
            for feature_idx in range(min(len(roi_features1), len(roi_features2))):
                features_list.append(
                    [roi_features1[feature_idx], roi_features2[feature_idx]]
                )
        features_list = np.array(features_list)
        new_shift = estimation_func(features_list, diagnostics=diagnostics)
        new_shift = np.round(new_shift).astype(np.int64)
        # if np.linalg.norm(new_shift, shift) <= tolerance:
        #     break
        shift = new_shift
    diagnostics["features1"] = hv.Points(features_list.swapaxes(0, 1)[0])
    diagnostics["features2"] = hv.Points(features_list.swapaxes(0, 1)[1])
    return shift
