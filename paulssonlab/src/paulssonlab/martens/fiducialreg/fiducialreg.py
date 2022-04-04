#!/usr/bin/env python
# -*- coding: utf-8 -*-
# fiducialreg.py

"""Generate transformation matrices from arrays of fiducial markers for
image registration.

:Author:
  `Talley Lambert <http://www.talleylambert.com>`_

:Organization:
  Cell Biology Microscopy Facility, Harvard Medical School

Acknowledgements
----------------
*   David Baddeley for gaussian fitting from python-microscopy
*   Andriy Myronenko for CPD algorithm
*   Siavash Khallaghi for python CPD implementation: pycpd

References
----------
(1) Point-Set Registration: Coherent Point Drift.  Andriy Mryonenko & Xubo Song
    https://arxiv.org/abs/0905.2635

Examples
--------
>>> from fiducialreg import FiducialCloud
>>> arrays = [arr1, arr2, arr3].  # three channels of fiducial marker stacks
>>> R = FiducialCloud(arrays, labels=[488, 560, 640])
>>> 560_to_488_rigid = R.get_tform_by_label(560, 488, mode='rigid')
>>> print(560_to_488_rigid)

# 560_to_488_rigid can then be used to transform a stack from 560 channel
# to a stack from the 488 channel...
>>> out = affine(im560, 560_to_488_rigid)
"""
from __future__ import division

from scipy import ndimage, stats
import numpy as np
import logging

logger = logging.getLogger(__name__)

logger.setLevel("DEBUG")
np.seterr(divide="ignore", invalid="ignore")


class RegistrationError(Exception):
    """Base class for fiducialreg errors"""

    pass


# TODO: try seperable gaussian filter instead for speed
def log_filter(img, blurxysigma=1, blurzsigma=2.5, mask=None):
    # sigma that works for 2 or 3 dimensional img
    sigma = [blurzsigma, blurxysigma, blurxysigma][-img.ndim :]
    # LOG filter image
    filtered_img = -ndimage.gaussian_laplace(img.astype("f"), sigma)
    # eliminate negative pixels
    filtered_img *= filtered_img > 0
    if mask is not None:
        filtered_img *= mask

    return filtered_img


def bead_centroids(img, labeled, nlabels):
    """get center of mass of each object"""
    return [ndimage.center_of_mass(img, labeled, l) for l in range(1, nlabels + 1)]


def get_thresh(im, mincount=None, steps=100):
    """intelligently find coordinates of local maxima in an image
    by searching a range of threshold parameters to find_local_maxima

    Accepts: variable number of input 2D arrays

    Returns:
    a tuple of sets of tuples ({(x,y),..},{(x,y),..},..) corresponding to
    local maxima in each image provided.  If nimages in == 1, returns a set

    """
    if im.ndim == 3:
        im = im.max(0)
    if mincount is None:
        mincount = 20
    threshrange = np.linspace(im.min(), im.max(), steps)
    object_count = [ndimage.label(im > t)[1] for t in threshrange]
    object_count = np.array(object_count)
    if mincount > object_count.max():
        raise RegistrationError(
            "Could not detect minimum number of beads specified ({}), found: {}".format(
                mincount, object_count.max()
            )
        )
    modecount = stats.mode(object_count[(object_count >= mincount)], axis=None).mode[0]
    logging.debug(
        "Threshold detected: {}".format(
            threshrange[np.argmax(object_count == modecount)]
        )
    )

    return threshrange[np.argmax(object_count == modecount)], modecount


def mad(arr, axis=None, method="median"):
    """Median/Mean Absolute Deviation: a "Robust" version of standard deviation.
    Indices variabililty of the sample.
    https://en.wikipedia.org/wiki/Median_absolute_deviation
    """
    if method == "median":
        return np.median(np.abs(arr - np.median(arr, axis)), axis)
    elif method == "mean":
        return np.mean(np.abs(arr - np.mean(arr, axis)), axis)
    else:
        raise ValueError("Unrecognized option for method: {}".format(method))


def get_closest_points(pc1, pc2):
    """returns the distance and index of the closest matching point in pc2
    for each point in pc1.

    len(nn) == len(pc1)

    can be used to eliminate points in pc2 that don't have a partner in pc1
    """
    pc1 = pc1.T
    pc2 = pc2.T
    d = [((pc2 - point) ** 2).sum(axis=1) for point in pc1]
    nn = [(np.min(p), np.argmin(p)) for p in d]
    return nn


def get_matching_points(pc1, pc2, method=None):
    """return modified point clouds such that every point in pc1 has a
    neighbor in pc2 that is within distance maxd
    """
    pc2neighbor_for_pc1 = np.array(get_closest_points(pc1, pc2))
    if method == "mean":
        mdist = np.mean(pc2neighbor_for_pc1, 0)[0]
        mdev = mad(pc2neighbor_for_pc1, 0, method="mean")[0]
    else:
        mdist = np.median(pc2neighbor_for_pc1, 0)[0]
        mdev = mad(pc2neighbor_for_pc1, 0, method="median")[0]
    passing = abs(pc2neighbor_for_pc1[:, 0] - mdist) < mdev * 4
    goodpc1 = pc1.T[passing]
    goodpc2 = pc2.T[pc2neighbor_for_pc1[:, 1][passing].astype("int")]

    return goodpc1.T, goodpc2.T


def mat2to3(mat2):
    """2D to 3D matrix:
    | a b c |       | a b 0 c |
    | d e f |  =>   | d e 0 f |
    | g h i |       | 0 0 1 0 |
                    | g h 0 i |
    """
    mat3 = np.eye(4)
    mat3[0:2, 0:2] = mat2[0:2, 0:2]
    mat3[3, 0:2] = mat2[2, 0:2]
    mat3[0:2, 3] = mat2[0:2, 2]
    mat3[3, 3] = mat2[2, 2]

    return mat3


# ### INFER TRANSFORMS ####


def infer_affine(X, Y, homo=1):
    """calculate affine transform which maps a set of points X onto Y

    X - 3xM XYZ points in starting coordinate system.
    Y - 3xM XYZ points in destination coordinate system.

    (The resulting transform will take points from X space and map
    them into Y space).
    """
    ndim = X.shape[0]
    if homo:
        X = np.vstack((X, np.ones((1, X.shape[1]))))

    # FIXME rcond=-1?
    affT = np.linalg.lstsq(X.T, Y.T, rcond=None)[0].T

    M = np.eye(ndim + 1)
    M[:ndim, :] = affT

    return M


def infer_rigid(X, Y, scale=False):
    n = X.shape[1]
    tVec = np.mean(Y - X, 1)

    # And the mean-corrected positions
    Mx = X - np.tile(np.mean(X, 1), (n, 1)).T
    My = Y - np.tile(np.mean(Y, 1), (n, 1)).T

    # Now solve for rotation matrix, using [1]
    CC = np.dot(My, Mx.T) / n
    U, _, V = np.linalg.svd(CC)
    F = np.eye(3)
    # Prevents reflection.
    F[2, 2] = np.linalg.det(np.dot(U, V))
    rMat = np.dot(np.dot(U, F), V)

    if scale:
        sigmaXsq = np.sum(Mx**2) / n
        scaling = np.trace(np.dot(rMat.T, CC)) / sigmaXsq
    else:
        scaling = 1
    # return rMat, tVec, rotCent, scaling
    # rotCent = np.mean(X, 1)

    # construct matrix
    ndim = X.shape[0]
    M = np.eye(ndim + 1)
    M[:ndim, :ndim] = rMat * scaling
    M[:ndim, ndim] = tVec

    return M


def infer_similarity(X, Y):
    return infer_rigid(X, Y, scale=True)


def infer_2step(X, Y):
    Yxyz = Y
    Yxy = Yxyz[:2]
    Xxy = X[:2]
    Xz = X[2:]
    T1 = infer_affine(Xxy, Yxy)
    M = mat2to3(T1)
    Xxy_reg = affineXF(Xxy, T1)
    Xxyz_reg = np.concatenate((Xxy_reg, Xz), axis=0)
    T2 = infer_similarity(Xxyz_reg, Yxyz)
    M[0:3, -1] += T2[0:3, -1]
    M[2, 2] *= T2[2, 2]

    return M


def infer_translation(X, Y):
    ndim = X.shape[0]
    M = np.eye(ndim + 1)
    M[0:ndim, -1] = np.mean(Y - X, 1)

    return M


# ### APPLY TRANSFORMS ####


def cart2hom(X):
    return np.vstack((X, np.ones((1, X.shape[1]))))


def intrinsicToWorld(intrinsicXYZ, dxy, dz):
    # def intrinsicToWorld(intrinsicXYZ, dxy, dz, worldStart=0.5):
    """where intrinsicXYZ is a 1x3 vector np.array([X, Y, Z])"""
    # NOTE by ATM: is the reasoning behind the 0.5 that it places
    # within the "center" of a pixel?
    # See: http://alvyray.com/Memos/CG/Microsoft/6_pixel.pdf
    # "A Pixel Is _Not_ a Little Square!"
    # If I understand Alvy Ray's argument, then the 0.5 nomenclature
    # is neither necessary nor desirable. It should be sufficient
    # to always assume a worldStart of 0.0, and to simply
    # multiply or divide by the dx & dz scaling factors.
    # By this reasoning, the point (0,0) remains at (0,0) regardless
    # of the underlying scaling, and other points are linearly
    # scaled. This means that intrinsic to world just multiplies
    # by dxy or dz, and world to intrinsic just divides by dxy or dz.
    # In fact, this is the procedure already invoked when fitting Gaussians!
    # Note that it might then be helpful to pre-store [dxy, dxy, dz]
    # as a numpy array within the class ~ would already be an ndarray,
    # and would not have to be created de novo every single invocation.
    if dxy == dz == 1:
        logger.warning("voxel size set at [1,1,1]... possibly unset")

    return intrinsicXYZ * np.array([dxy, dxy, dz])
    # return worldStart + (intrinsicXYZ - 0.5) * np.array([dxy, dxy, dz])


def worldToInstrinsic(worldXYZ, dxy, dz):
    # def worldToInstrinsic(worldXYZ, dxy, dz, worldStart=0.5):
    """where XYZ coord is a 1x3 vector np.array([X, Y, Z])"""
    if dxy == dz == 1:
        logger.warning("voxel size set at [1,1,1]... possibly unset")

    # NOTE see above about doing without the +/- 0.5 conventions
    return worldXYZ / np.array([dxy, dxy, dz])
    # return 0.5 + (worldXYZ - worldStart) / np.array([dxy, dxy, dz])


def affineXF(X, T, invert=False):
    ndim = X.shape[0]
    X = np.vstack((X, np.ones((1, X.shape[1]))))
    if not invert:
        return np.dot(T, X)[:ndim, :]
    else:
        return np.dot(np.linalg.inv(T), X)[:ndim, :]


def rigidXF(X, rMat, tVec, rotCent=None, scaling=1, invert=False):
    xlen = X.shape[1]

    if rotCent is None:
        rotCent = np.mean(X, 1)

    X - np.tile(rotCent, (xlen, 1)).T

    if not invert:
        Y = np.dot(rMat, X - np.tile(rotCent, (xlen, 1)).T)
        Y *= scaling
        Y += np.tile(rotCent + tVec, (xlen, 1)).T
    else:
        Y = np.dot(np.linalg.inv(rMat), X - np.tile(rotCent + tVec, (xlen, 1)).T)
        Y /= scaling
        Y += np.tile(rotCent, (xlen, 1)).T

    return Y


def translateXF(X, T, invert=False):
    T = np.tile(T[0:3, 3], (X.shape[1], 1))
    if not invert:
        return X + T.T
    else:
        return X - T.T
