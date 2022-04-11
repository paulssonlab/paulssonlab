#!/usr/bin/env python

from __future__ import division

from scipy import ndimage
from os import path as osp
import numpy as np
import logging
import json

from .fiducialreg import get_thresh, intrinsicToWorld
from .gaussfit import fit_gaussian

# using Qt5Agg causes "window focus loss" in interpreter for some reason
import matplotlib

# matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)
logger.setLevel("WARNING")


class lazyattr(object):
    """Lazy object attribute whose value is computed on first access."""

    __slots__ = ("func",)

    def __init__(self, func):
        self.func = func

    def __get__(self, instance, owner):
        if instance is None:
            return self
        value = self.func(instance)
        if value is NotImplemented:
            return getattr(super(owner, instance), self.func.__name__)
        setattr(instance, self.func.__name__, value)
        return value


class FiducialCloud(object):
    """Generate a 3D point cloud of XYZ locations of fiducial markers

    Points are extracted from an image of (e.g.) beads by 3D gaussian fitting.
    Creates an instance for a single wavelength.  Multiple Clouds can be related
    to each other using a :obj:`CloudSet` object.

    Args:
        data (:obj:`np.ndarray`, :obj:`str`): 3D numpy array or
            path to LLS registration folder experiment
        dz (:obj:`float`): Z-step size
        dx (:obj:`float`): Pixel size in XY
        blurxysig (:obj:`float`): XY sigma for pre-fitting blur step
        blurzsig (:obj:`float`): Z sigma for pre-fitting blur step
        threshold ('auto', :obj:`int`): intensity threshold when detecting beads.
            by default, it will autodetect threshold using mincount parameter.
        mincount (:obj:`int`): minimum number of beads expected in the image,
            used when autodetecting threshold.
        imref (:obj:`fiducialreg.imref.imref3d`): spatial referencing object.
        filtertype ({'blur', 'log'}): type of blur to use prior to bead detection.
            log = laplacian of gaussian
            blur = gaussian

    """

    def __init__(
        self,
        data=None,
        dz=1,
        dx=1,
        blurxysig=1,
        blurzsig=2.5,
        threshold=None,
        mincount=None,
        imref=None,
        filtertype="blur",
    ):
        # data is a numpy array or filename
        self.data = None
        if data is not None:
            if isinstance(data, str) and osp.isfile(data):
                # fc = FiducialCloud('/path/to/file')
                try:
                    import tifffile as tf

                    self.filename = osp.basename(data)
                    self.data = tf.imread(data).astype("f")
                except ImportError:
                    raise ImportError(
                        "The tifffile package is required to read a "
                        "filepath into an array."
                    )
            elif isinstance(data, np.ndarray):
                # fc = FiducialCloud(np.ndarray)
                self.data = data
            else:
                raise ValueError(
                    "Input to Registration must either be a "
                    "filepath or a numpy arrays"
                )
        self.dx = dx
        self.dz = dz
        self.blurxysig = blurxysig
        self.blurzsig = blurzsig
        self.threshold = threshold if threshold is not None else "auto"
        self._mincount = mincount
        self.imref = imref
        self.coords = None
        self.filtertype = filtertype

        logger.debug("New fiducial cloud created with dx: {},  dz: {}".format(dx, dz))
        if self.data is not None:
            self.update_coords()

    def has_data(self):
        return isinstance(self.data, np.ndarray)

    @property
    def mincount(self):
        return self._mincount

    @mincount.setter
    def mincount(self, value):
        self._mincount = value
        self.update_coords()
        logger.info("found {} spots".format(self.count))

    @property
    def count(self):
        if self.coords is not None and self.coords.shape[1] > 0:
            return self.coords.shape[1]
        else:
            return 0

    @lazyattr
    def filtered(self):
        if self.data is not None:
            if self.filtertype == "log":
                return log_filter(self.data, self.blurxysig, self.blurzsig)
            else:
                ### EDITED BY ATM
                # This might have been included because it was packaged with
                # LLSPY, which uses other GPU offloading?
                # For our use case definitely not worth the hassle!
                # from gputools import blur
                # sigs = np.array([self.blurzsig, self.blurxysig, self.blurxysig]) * 2
                # import warnings
                # with warnings.catch_warnings():
                # warnings.simplefilter("ignore")
                # out = blur(self.data, sigs)
                # return out

                # NOTE is doing 2D gaussian blurs
                # Seems to work, unclear if results are identical
                # from gputools.
                # NOTE this will fail if only have 2D images??
                # in which case, zap the blurzsig bit
                from scipy.ndimage import gaussian_filter as blur

                sigs = np.array([self.blurzsig, self.blurxysig, self.blurxysig]) * 2
                out = blur(self.data, sigs)

                return out
        else:
            return None

    def autothresh(self, mincount=None):
        if mincount is None:
            mincount = self._mincount
        return get_thresh(self.filtered, mincount)[0]

    def update_coords(self, thresh=None):
        if self.filtered is None:
            return
        if thresh is None:
            thresh = self.threshold
        if thresh == "auto" or not thresh:
            thresh = self.autothresh()
        try:
            thresh = float(thresh)
            assert thresh > 0
        except Exception:
            raise RegistrationError(
                "Threshold must be number greater than 0.  got: {}".format(thresh)
            )
        logger.debug("Update_coords using threshold: {}".format(thresh))

        # Assign unique labels
        labeled = ndimage.label(self.filtered > thresh)[0]

        # Find objects
        objects = ndimage.find_objects(labeled)

        # Append new spots to this list
        coords_xyz = []

        # To normalize spot coordinates
        denominator = np.array([self.dx, self.dx, self.dz])

        # Iterate the chunks, fit gaussians
        for chunknum, chunk in enumerate(objects):
            try:
                # Get the xyz coordinates of this spot by fitting to a 3D gaussian
                # FIXME: pass sigmas to wx and wz parameters
                # TODO: filter by bead intensity as well to reject bright clumps

                # F: 0: x, 1: y, 2: z
                F = fit_gaussian(data=self.data, dx=self.dx, dz=self.dz, key=chunk)

                # Convert from "world" space back into "intrinsic" space
                # NOTE how come here we don't need to do the whole intrinsic -> world
                # conversion, which includes a 0.5 offset & a WorldStart value?
                # NOTE previous code had "realspace=True" as a default value, but actually called
                # with realspace=False. If need be, we could optionally not normalize.
                F = np.divide(F, denominator)

                # Verify that the resulting fit is within the bounds of the image
                if (
                    (F[0] < self.data.shape[2])
                    and (F[0] > 0)
                    and (F[1] < self.data.shape[1])
                    and (F[1] > 0)
                    and (F[2] < self.data.shape[0])
                    and (F[2] > 0)
                ):
                    logger.debug("Added spot number {}, at {}".format(chunknum, chunk))
                    coords_xyz.append(F)

            except Exception:
                logging.warning("Skipped spot number {}, at {}".format(chunknum, chunk))

        # Transpose & store into the cloud object
        self.coords = np.array(coords_xyz).T

        if not len(self.coords):
            logging.warning(
                "PointCloud has no points! {}".format(
                    self.filename if "filename" in dir(self) else ""
                )
            )
        else:
            logger.info("Registration found {} objects".format(self.count))

    @property
    def coords_inworld(self):
        """return coordinate list in world coordinates based on voxel size"""
        return intrinsicToWorld(self.coords.T, self.dx, self.dz).T

    def show(self, withimage=True, filtered=True):
        if withimage and self.filtered is not None:
            if filtered:
                im = self.filtered.max(0)
            else:
                im = self.data.max(0)
            plt.imshow(im, cmap="gray", vmax=im.max() * 0.7)
        if self.count:
            plt.scatter(self.coords[0], self.coords[1], c="red", s=5)

    def toJSON(self):
        D = self.__dict__.copy()
        D.pop("filtered", None)
        D.pop("data", None)
        D["coords"] = self.coords.tolist()
        return json.dumps(D)

    def fromJSON(self, Jstring):
        logger.debug("Recovering Fidicuals from JSON string")
        J = json.loads(Jstring)
        for k, v in J.items():
            if not k == "coords":
                logger.debug("Setting attribute {} to {}".format(k, v))
            else:
                logger.debug("Populating coordinates from JSON string")
            setattr(self, k, v)
        self.coords = np.array(self.coords)
        return self
