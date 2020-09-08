import numpy as np
from scipy.optimize import fmin


def get_poly_coeffs(g):
    npts = np.size(g)
    x = np.linspace(1, -1, npts)
    H = np.zeros((npts, 2))
    H[:, 0] = 1
    H[:, 1] = x
    coeffs = np.dot(np.linalg.pinv(H), g)
    remainder = np.sum((g - coeffs[0] - (coeffs[1] * x)) ** 2)
    return (coeffs[0], coeffs[1], remainder)


def get_poly_coeffs2(g):
    npts = np.size(g)
    x = -1.0 + np.arange(npts) * 2.0 / (npts - 1.0)
    H = np.zeros((np.size(g), 3))
    H[:, 0] = 1.0
    H[:, 1] = x
    H[:, 2] = x**2
    coeffs = np.dot(np.linalg.pinv(H), g)
    remainder = np.sum((g - coeffs[0] - (coeffs[1] * x) - (coeffs[2] * x**2)) ** 2)
    return (coeffs[0], coeffs[1], coeffs[2], remainder)


def nonlinearity(x):
    g = np.gradient(np.gradient(x))
    return np.sqrt(np.sum(g**2))


def nonlinearity_error(x):
    grad = np.gradient(x)
    if np.max(grad) * np.min(grad) < 0.0:
        grad = grad - np.min(grad)
    g = np.gradient(grad)
    NL_err = np.sum(g**2)
    return NL_err


def spectral_sampling_ratio(
    w, delta_spectrum, shift=False, sampling_domain="Wavelength"
):
    if sampling_domain == "wavenumber":
        w = 1000 / w
    grad = np.gradient(delta_spectrum) / np.gradient(w)
    maxval = np.max(grad)
    minval = np.min(grad)

    ## If both maxval and minval are negative, then we should swap the ratio to
    ## get ssr > 1.
    if np.abs(maxval) < np.abs(minval):
        ssr = minval / maxval
    else:
        ssr = maxval / minval

    ## For the optimizer, we don't want to allow the gradient to cross zero, so
    ## we can penalize the ssr value by shifting the curve up just enough to
    ## prevent zero-crossing.
    if shift and (ssr < 0.0):
        grad = grad - np.min(grad)
        if np.abs(maxval) < np.abs(minval):
            ssr = minval / maxval
        else:
            ssr = maxval / minval

    return ssr


def beam_compression(thetas):
    center_wavelength_idx = thetas.shape[1] // 2
    abscostheta = np.abs(np.cos(thetas))
    K = np.product(abscostheta[1::2], axis=0) / np.product(abscostheta[:-1:2], axis=0)
    return K[center_wavelength_idx]


def calc_delta_singlet(angles, n1):
    theta0 = 0.0
    beta = angles[0]
    gamma = angles[1]
    alpha = gamma - beta
    theta1 = theta0 - beta
    theta1p = np.arcsin((1 / n1) * np.sin(theta1))
    theta2 = theta1p - alpha
    theta2p = np.arcsin(n1 * np.sin(theta2))
    theta3 = theta2p + gamma
    delta_spectrum = theta0 - theta3

    thetas = np.array([theta1 * np.ones(np.size(n1)), theta1p, theta2, theta2p, theta3])

    return (delta_spectrum, thetas)


def err_singlet(angles, n1, deltaT_target, angle_limit):
    (delta_spectrum, thetas) = calc_delta_singlet(angles, n1)
    deltaT_err = ((delta_spectrum.max() - delta_spectrum.min()) - deltaT_target) ** 2
    merit_err = deltaT_err

    ## If TIR occurs in the design (producing NaNs in the spectrum), then give a
    ## hard error: return a large error which has nothing to do with the (invalid)
    ## performance data.
    if np.any(np.isnan(delta_spectrum)):
        return 0.1 * np.sum(np.isnan(delta_spectrum))
    if np.any(np.iscomplex(delta_spectrum)):
        return 0.1 * np.sum(np.iscomplex(delta_spectrum))

    ## If any of the angles of incidence exceed the limit, increase the error value.
    if np.any(np.abs(np.array([thetas[0], thetas[2]])) > angle_limit):
        num_bad = np.sum(np.abs(np.array([thetas[0], thetas[2]])) > angle_limit)
        merit_err += num_bad / np.size(delta_spectrum)

    return merit_err


def calc_delta_doublet(angles, n1, n2, system="doublet"):
    theta0 = 0.0
    alpha1 = angles[0]
    alpha2 = angles[1]
    beta1 = alpha1 + (0.5 * alpha2)
    theta1 = theta0 + beta1
    theta1p = np.arcsin((1.0 / n1) * np.sin(theta1))
    theta2 = theta1p - alpha1
    theta2p = np.arcsin((n1 / n2) * np.sin(theta2))
    theta3 = theta2p - alpha2
    if system == "doublet":
        theta3p = np.arcsin(n2 * np.sin(theta3))
        theta4 = theta3p + 0.5 * alpha2
        delta_spectrum = theta0 - theta4
        thetas = np.array(
            [
                theta1 * np.ones(np.size(n1)),
                theta1p,
                theta2,
                theta2p,
                theta3,
                theta3p,
                theta4,
            ]
        )
    elif system == "double-amici":
        theta3p = np.arcsin((n2 / n1) * np.sin(theta3))
        theta4 = theta3p - alpha1
        theta4p = np.arcsin(n1 * np.sin(theta4))
        theta5 = theta4p + beta1
        delta_spectrum = theta0 - theta5
        thetas = np.array(
            [
                theta1 * np.ones(np.size(n1)),
                theta1p,
                theta2,
                theta2p,
                theta3,
                theta3p,
                theta4,
                theta4p,
                theta5,
            ]
        )
    else:
        raise NotImplementedError
    return (delta_spectrum, thetas)


def err_doublet(
    angles, n1, n2, deltaC_target, deltaT_target, merit, angle_limit, system="doublet"
):
    (delta_spectrum, thetas) = calc_delta_doublet(angles, n1, n2, system=system)

    ## If TIR occurs in the design (producing NaNs in the spectrum), then give a
    ## hard error: return a large error which has nothing to do with the (invalid)
    ## performance data.
    if np.any(np.isnan(delta_spectrum)):
        return 0.1 * np.sum(np.isnan(delta_spectrum))

    if merit == "chromaticity":
        deltaC_err = (delta_spectrum[np.size(delta_spectrum) // 2] - deltaC_target) ** 2
        chromat = delta_spectrum.max() - delta_spectrum.min()
        merit_err = deltaC_err + chromat
    elif merit == "dev_linearity":
        NL_err = 100.0 * nonlinearity_error(delta_spectrum)
        deltaT_err = ((delta_spectrum[0] - delta_spectrum[-1]) - deltaT_target) ** 2
        merit_err = deltaT_err + NL_err
    elif merit == "dev_linearityK":
        NL_err = 100.0 * nonlinearity_error(delta_spectrum)
        deltaT_err = ((delta_spectrum[0] - delta_spectrum[-1]) - deltaT_target) ** 2
        K = beam_compression(thetas)
        K_err = 0.01 * (K - 1.0) ** 2
        merit_err = deltaT_err + NL_err + K_err
    else:
        deltaC_err = (delta_spectrum[np.size(delta_spectrum) // 2] - deltaC_target) ** 2
        deltaT_err = ((delta_spectrum[0] - delta_spectrum[-1]) - deltaT_target) ** 2
        merit_err = deltaC_err + deltaT_err

    ## If any of the angles of incidence exceed the limit, increase the error value.
    if np.any(np.abs(np.array([thetas[0], thetas[2], thetas[4]])) > angle_limit):
        num_bad = np.sum(
            np.abs(np.array([thetas[0], thetas[2], thetas[4]])) > angle_limit
        )
        merit_err += num_bad / np.size(delta_spectrum)

    return merit_err


def calc_delta_triplet(angles, n1, n2, n3):
    theta0 = 0.0
    lim = pi / 2.1
    alpha1 = angles[0]
    alpha2 = angles[1]
    alpha3 = angles[2]
    beta1 = alpha1 + (0.5 * alpha2)
    theta1 = theta0 + beta1
    theta1p = arcsin((1.0 / n1) * sin(theta1))
    theta2 = theta1p - alpha1
    theta2p = arcsin((n1 / n2) * sin(theta2))
    theta3 = theta2p - alpha2
    theta3p = arcsin((n2 / n3) * sin(theta3))
    theta4 = theta3p - alpha3
    theta4p = arcsin(n3 * sin(theta4))
    gamma3 = alpha3 + (0.5 * alpha2)
    theta5 = theta4p + gamma3
    delta_spectrum = theta0 - theta5

    thetas = array(
        [
            theta1 * ones(size(n1)),
            theta1p,
            theta2,
            theta2p,
            theta3,
            theta3p,
            theta4,
            theta4p,
            theta5,
        ]
    )

    return (delta_spectrum, thetas)


def err_damici(angles, n1, n2, n3, deltaC_target, deltaT_target, merit, angle_limit):
    (delta_spectrum, thetas) = calc_delta_triplet(angles, n1, n2, n3)

    merit_err = 0.0

    ## If TIR occurs in the design (producing NaNs in the spectrum), then give a
    ## hard error: return a large error which has nothing to do with the (invalid)
    ## performance data.
    if np.any(np.isnan(delta_spectrum)):
        return 0.1 * np.sum(np.isnan(delta_spectrum))

    ## If any of the angles of incidence exceed the limit, increase the error value.
    if np.any(
        np.abs(np.array([thetas[0], thetas[2], thetas[4], thetas[6]])) > angle_limit
    ):
        num_bad = np.sum(
            np.abs(np.array([thetas[0], thetas[2], thetas[4], thetas[6]])) > angle_limit
        )
        merit_err += num_bad / np.size(delta_spectrum)

    if merit == "linearity":
        NL_err = 100.0 * nonlinearity_error(delta_spectrum)
        deltaC_err = (delta_spectrum[size(delta_spectrum) // 2] - deltaC_target) ** 2
        deltaT_err = (
            (delta_spectrum.max() - delta_spectrum.min()) - deltaT_target
        ) ** 2
        merit_err += deltaC_err + deltaT_err + NL_err
    elif merit == "linearityK":
        NL_err = 25.0 * nonlinearity_error(delta_spectrum)
        deltaC_err = (delta_spectrum[size(delta_spectrum) // 2] - deltaC_target) ** 2
        deltaT_err = (
            (delta_spectrum.max() - delta_spectrum.min()) - deltaT_target
        ) ** 2
        K = beam_compression_triplet(thetas)
        K_err = 0.0025 * (K - 1.0) ** 2
        merit_err += deltaC_err + deltaT_err + NL_err + K_err
    elif merit == "dev_linearity":
        NL_err = 25.0 * nonlinearity_error(delta_spectrum)
        deltaT_err = (
            (delta_spectrum.max() - delta_spectrum.min()) - deltaT_target
        ) ** 2
        merit_err += deltaT_err + NL_err
    elif merit == "dev_linearityK":
        NL_err = 25.0 * nonlinearity_error(delta_spectrum)
        deltaT_err = (
            (delta_spectrum.max() - delta_spectrum.min()) - deltaT_target
        ) ** 2
        K = beam_compression_triplet(thetas)
        K_err = 0.0025 * (K - 1.0) ** 2
        merit_err += deltaT_err + NL_err + K_err
    elif merit == "old_linearity":
        (mean_delta, delta1, remainder) = get_poly_coeffs(delta_spectrum)
        NL_err = remainder / np.size(delta_spectrum)
        deltaC_err = (delta_spectrum[np.size(delta_spectrum) // 2] - deltaC_target) ** 2
        deltaT_err = (
            (delta_spectrum.max() - delta_spectrum.min()) - deltaT_target
        ) ** 2
        merit_err += deltaC_err + deltaT_err + NL_err
    elif merit == "chromaticity":
        deltaC_err = (delta_spectrum[np.size(delta_spectrum) // 2] - deltaC_target) ** 2
        chromat = delta_spectrum.max() - delta_spectrum.min()
        merit_err += deltaC_err + chromat
    elif merit == "angle_sum":
        anglesum = 0.001 * deltaT_target * sum(angles**2)
        deltaC_err = (delta_spectrum[np.size(delta_spectrum) // 2] - deltaC_target) ** 2
        deltaT_err = (
            (delta_spectrum.max() - delta_spectrum.min()) - deltaT_target
        ) ** 2
        merit_err += deltaC_err + deltaT_err + anglesum
    elif merit == "second-order":
        (mean_delta, delta1, delta2, remainder) = get_poly_coeffs2(delta_spectrum)
        secondorder = 0.001 / delta2**2
        deltaC_err = (delta_spectrum[np.size(delta_spectrum) // 2] - deltaC_target) ** 2
        deltaT_err = (
            (delta_spectrum.max() - delta_spectrum.min()) - deltaT_target
        ) ** 2
        merit_err += deltaC_err + deltaT_err + secondorder

    if np.any(np.abs(angles) < 1.0 * pi / 180.0):
        too_thin = np.abs(angles) - 1.0
        too_thin[(np.where(np.abs(angles) > 1.0 * pi / 180.0))[0]] = 0.0
        merit_err += 0.25 * np.sum(too_thin**2)

    return merit_err


def calc_delta_triplet(angles, n1, n2):
    theta0 = 0.0
    alpha1 = angles[0]
    alpha2 = angles[1]
    beta1 = alpha1 + (0.5 * alpha2)
    theta1 = theta0 + beta1
    theta1p = np.arcsin((1.0 / n1) * np.sin(theta1))
    theta2 = theta1p - alpha1
    theta2p = np.arcsin((n1 / n2) * np.sin(theta2))
    theta3 = theta2p - alpha2
    theta3p = np.arcsin((n2 / n1) * np.sin(theta3))
    theta4 = theta3p - alpha1
    theta4p = np.arcsin(n1 * np.sin(theta4))
    theta5 = theta4p + beta1
    delta_spectrum = theta0 - theta5
    thetas = np.array(
        [
            theta1 * np.ones(np.size(n1)),
            theta1p,
            theta2,
            theta2p,
            theta3,
            theta3p,
            theta4,
            theta4p,
            theta5,
        ]
    )
    return (delta_spectrum, thetas)


def err_triplet(angles, n1, n2, deltaC_target, deltaT_target, merit, angle_limit):
    (delta_spectrum, thetas) = calc_delta_triplet(angles, n1, n2)

    merit_err = 0.0

    ## If TIR occurs in the design (producing NaNs in the spectrum), then give a
    ## hard error: return a large error which has nothing to do with the (invalid)
    ## performance data.
    if any(isnan(delta_spectrum)):
        return 0.1 * sum(isnan(delta_spectrum))

    if merit == "chromaticity":
        deltaC_err = (delta_spectrum[size(delta_spectrum) // 2] - deltaC_target) ** 2
        chromat = delta_spectrum.max() - delta_spectrum.min()
        merit_err = deltaC_err + chromat
    elif merit == "dev_linearity":
        NL_err = 100.0 * nonlinearity_error(delta_spectrum)
        deltaT_err = ((delta_spectrum[0] - delta_spectrum[-1]) - deltaT_target) ** 2
        merit_err = deltaT_err + NL_err
    elif merit == "dev_linearityK":
        NL_err = 100.0 * nonlinearity_error(delta_spectrum)
        deltaT_err = ((delta_spectrum[0] - delta_spectrum[-1]) - deltaT_target) ** 2
        K = beam_compression(thetas)
        K_err = 0.01 * (K - 1.0) ** 2
        merit_err = deltaT_err + NL_err + K_err
    else:
        deltaC_err = (delta_spectrum[size(delta_spectrum) // 2] - deltaC_target) ** 2
        deltaT_err = ((delta_spectrum[0] - delta_spectrum[-1]) - deltaT_target) ** 2
        merit_err = deltaC_err + deltaT_err

    ## If any of the angles of incidence exceed the limit, increase the error value.
    if any(abs(array([thetas[0], thetas[2], thetas[4]])) > angle_limit):
        num_bad = sum(abs(array([thetas[0], thetas[2], thetas[4]])) > angle_limit)
        merit_err += num_bad / size(delta_spectrum)

    return merit_err


def optimize_singlet(w, n1, dispersion_target, sampling_domain="wavelength"):
    angle_limit = np.deg2rad(85)
    # dispersion_target = 1  ## in degrees for the full spectral width
    deltaT_target = dispersion_target  ## in radians per full spectral width
    # w, n1 = glasscat
    initial_angles = np.deg2rad([-5, 5])
    results = fmin(
        err_singlet, initial_angles, args=(n1, deltaT_target, angle_limit), disp=0
    )
    angles = np.asarray([results[0], results[1]])
    (delta_spectrum, thetas) = calc_delta_singlet(angles, n1)
    if np.any(np.isnan(delta_spectrum)):
        raise ValueError

    (deltaM, delta1, remainder) = get_poly_coeffs(delta_spectrum)
    deltaC = delta_spectrum[len(delta_spectrum) // 2]
    deltaT = delta_spectrum.max() - delta_spectrum.min()
    beta = angles[0]
    gamma = angles[1]
    alpha = gamma - beta
    nonlin = np.sqrt(remainder) / abs(delta1)
    NL = 10000.0 * nonlinearity(delta_spectrum)
    SSR = spectral_sampling_ratio(w, delta_spectrum)
    return (beta, gamma, alpha, deltaC, deltaM, delta1, deltaT, nonlin, NL, SSR)


def optimize_doublet(
    w,
    n1,
    n2,
    dispersion_target,
    system="doublet",
    merit="dev_linearityK",
    sampling_domain="wavelength",
):
    if system not in ("doublet", "double-amici"):
        raise ValueError(f"unknown system: {system}")
    if merit == "dev_linearity" or merit == "dev_linearityK":
        deviation_target = np.nan
    angle_limit = np.deg2rad(85)
    deltaC_target = deviation_target
    deltaT_target = dispersion_target
    if not np.isnan(deviation_target):
        # initial_angles = array([(deltaC_target/4.0)+deltaT_target, (deltaC_target/4.0)-deltaT_target])
        initial_angles = np.array([1, 1]) * deltaC_target / 4.0
        initial_angles += np.array([1, -1]) * deltaT_target
    else:
        # initial_angles = array([deltaT_target, -deltaT_target])
        initial_angles = np.array([1, -1]) * deltaT_target
    results = fmin(
        err_doublet,
        initial_angles,
        args=(n1, n2, deltaC_target, deltaT_target, merit, angle_limit, system),
        disp=0,
    )
    angles = np.asarray([results[0], results[1]])
    (delta_spectrum, thetas) = calc_delta_doublet(angles, n1, n2, system=system)
    if np.any(np.isnan(thetas)):
        raise ValueError

    (deltaM, delta1, delta2, remainder) = get_poly_coeffs2(delta_spectrum)
    if sampling_domain == "wavelength":
        delta1 = -delta1
    if delta1 < 0.0:  ## flip over the prism if the dispersion is inverted
        angles = -1 * angles
        delta1 = -1 * delta1
    deltaM = np.mean(delta_spectrum)
    deltaC = delta_spectrum[len(delta_spectrum) // 2]
    deltaT = delta_spectrum.max() - delta_spectrum.min()
    alpha1 = angles[0]
    alpha2 = angles[1]
    NL = 10000.0 * nonlinearity(delta_spectrum)
    SSR = spectral_sampling_ratio(w, delta_spectrum, sampling_domain=sampling_domain)
    K = beam_compression(thetas)
    ## Now use the two-coeff version of the polynomial fit; the nonlin "remainder"
    ## should not subtract out the second-order coeff, as the three-coeff version does.
    (temp1, temp2, remainder) = get_poly_coeffs(delta_spectrum)
    dw = np.abs(np.mean(np.gradient(w)))
    nonlin = np.sqrt(remainder) * dw
    chromat = 100.0 * (delta_spectrum.max() - delta_spectrum.min())
    return (
        alpha1,
        alpha2,
        deltaC,
        deltaT,
        NL,
        SSR,
        deltaM,
        K,
        delta1 * 2,
        delta2 * 2,
        nonlin,
        chromat,
    )
