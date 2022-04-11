#!/usr/bin/env python

from __future__ import division

from scipy import optimize
import numpy as np


def calculate_sigma(electrons_per_ADU, TrueEMGain, NoiseFactor, ReadNoise, dataROI):
    """
    should use gain and noise map from camera Parameters
    for now assume uniform noise characteristcs and sCMOS read noise
    """
    # estimate noise as read noise plus poisson noise
    sigma = (
        np.sqrt(
            ReadNoise**2
            + NoiseFactor**2 * electrons_per_ADU * TrueEMGain * np.maximum(dataROI, 1)
        )
        / electrons_per_ADU
    )

    return sigma


def f_Gauss3d(p, X, Y, Z):
    """3D PSF model function with constant background
    parameter vector [A, x0, y0, z0, background]
    """
    A, x0, y0, z0, wxy, wz, b = p
    # return A*scipy.exp(-((X-x0)**2 + (Y - y0)**2)/(2*s**2)) + b

    return (
        A
        * np.exp(
            -((X - x0) ** 2 + (Y - y0) ** 2) / (2 * wxy**2)
            - ((Z - z0) ** 2) / (2 * wz**2)
        )
        + b
    )


def weightedMissfitF(p, fcn, data, weights, *args):
    """Helper function which evaluates a model function (fcn) with parameters (p)
    and additional arguments (*args) and compares this with measured data (data),
    scaling with precomputed weights corresponding to the errors in the measured
    data (weights).
    """
    model = fcn(p, *args)
    model = model.ravel()

    return (data - model) * weights


def FitModelWeighted(modelFcn, startParameters, data, sigmas, *args):
    """
    NOTE: full_output is useful for calculating the fit errors.
    infodict['fvec']
    """
    return optimize.leastsq(
        weightedMissfitF,
        startParameters,
        (modelFcn, data.ravel(), (1.0 / sigmas).astype("f").ravel()) + args,
        full_output=True,
    )


# Edited by ATM. Trying to make the code less complex by eliminating a class,
# replacing it with just a function (for running the fits),
# and another class, replacing it with a dict (for returning the results).
# The main functionality that the result seemed to provide was a mechanism
# for automatically normalizing by dx & dz, but that hardly seems worth
# the added complexity!


def fit_gaussian(data, dx, dz, key, wx=0.17, wz=0.37):
    """return gaussian fit of a 3D roi defined by a 3-tuple of slices"""
    zslice, yslice, xslice = key
    # cut region out of data stack
    dataROI = data[zslice, yslice, xslice].astype("f")

    # generate grid to evaluate function on
    Z, Y, X = np.mgrid[zslice, yslice, xslice]

    # adjust for voxel size
    # NOTE how come here we don't need to do the whole intrinsic -> world
    # conversion, which includes a 0.5 offset & a WorldStart value?
    X = dx * X
    Y = dx * Y
    Z = dz * Z

    # amplitude
    A = dataROI.max() - dataROI.min()

    # subtract background
    drc = dataROI - dataROI.min()

    drc = np.maximum(drc - drc.max() / 2, 0)

    # normalize sum to 1
    drc = drc / drc.sum()

    x0 = (X * drc).sum()
    y0 = (Y * drc).sum()
    z0 = (Z * drc).sum()

    startParameters = [3 * A, x0, y0, z0, wx, wz, dataROI.min()]

    # NOTE Does this estimate the amount of noise from the camera?
    sigma = calculate_sigma(
        electrons_per_ADU=0.5,
        TrueEMGain=1,
        NoiseFactor=1,
        ReadNoise=1.2,
        dataROI=dataROI,
    )

    # Calculate the fit & some other information:
    # covariance matrix, infodict, mesgl, resCode
    (res1, cov_x, infodict, mesg1, resCode) = FitModelWeighted(
        f_Gauss3d, startParameters, dataROI, sigma, X, Y, Z
    )
    # misfit = (infodict['fvec']**2).sum()  # nfev is the number of function calls

    # Also calculate the "errors" of the fit
    # NOTE the errors are never used! So let's comment this out for now.
    # fitErrors = None
    # try:
    # fitErrors = np.sqrt(
    # np.diag(cov_x)
    # * (infodict["fvec"] * infodict["fvec"]).sum()
    # / (dataROI.size - len(res1))
    # )
    # except Exception:
    # pass

    # Return an ndarray of size 3, with xyz values
    return_val = np.array([res1[1], res1[2], res1[3]])

    # NOTE we could also return a dict.
    # Note that fiducialcloud only really cares about x, y, & z.
    # return_val = {
    # "A"          : res1[0],
    # "x"          : res1[1],
    # "y"          : res1[2],
    # "z"          : res1[3],
    # "wxy"        : res1[4],
    # "wz"         : res1[5],
    # "background" : res1[6],
    # "dx"         : dx,
    # "dz"         : dz,
    # "key"        : key,
    # "resCode"    : resCode,
    # "fitErrors"  : fitErrors
    # }

    return return_val
