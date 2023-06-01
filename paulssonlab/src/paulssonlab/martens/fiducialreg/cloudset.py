#!/usr/bin/env python

from __future__ import division

import itertools
import json
import logging
from os import path as osp

# using Qt5Agg causes "window focus loss" in interpreter for some reason
import matplotlib

# matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas

from .cpd import CPDaffine, CPDrigid, CPDsimilarity, cpd_2step
from .fiducialcloud import FiducialCloud
from .fiducialreg import (
    get_matching_points,
    infer_2step,
    infer_affine,
    infer_rigid,
    infer_similarity,
    infer_translation,
)

logger = logging.getLogger(__name__)
logger.setLevel("DEBUG")


class CloudSet(object):
    """Creates a set of fiducial clouds for the purpose of estimating transform

    The main method is the tform() method to generate a transformation matrix
    mapping one cloud in the set to another.

    Args:
        data (Iterable[:obj:`str`, :obj:`np.ndarray`]): list of 3D numpy arrays,
            or a list of filepath strings.  Strings will be handed to
            :obj:`FiducialCloud` and read using tifffile
        labels (Iterable[:obj:`str`], optional):  List of labels (e.g. wave) naming each
            fiducial cloud in the set.  This allows registration using labels
            instead of the index of the original dataset provided to :obj:`CloudSet`

        **kwargs: extra keyword arguments are passed to the :obj:`FiducialCloud`
            constructor.

    """

    def __init__(
        self,
        data=None,
        labels=None,
        dx=1,
        dz=1,
        mincount=None,
        threshold=None,
        **kwargs,
    ):
        self.dx = dx
        self.dz = dz
        if data is not None:
            if not isinstance(data, (list, tuple, set)):
                raise ValueError(
                    "CloudSet expects a list of np.ndarrays or " "filename strings"
                )
            if labels is not None:
                if len(labels) != len(data):
                    raise ValueError(
                        "Length of optional labels list must match "
                        "length of the data list"
                    )
                self.labels = labels
            else:
                self.labels = ["ch" + str(i) for i in range(len(data))]
            self.N = len(data)
            self.clouds = []
            for i, d in enumerate(data):
                logger.info(
                    "Creating FiducalCloud for label: {}".format(self.labels[i])
                )
                self.clouds.append(
                    FiducialCloud(
                        d,
                        dx=self.dx,
                        dz=self.dz,
                        threshold=threshold,
                        mincount=mincount,
                        **kwargs,
                    )
                )
        else:
            self.clouds = []
            self.N = 0

    def toJSON(self):
        return json.dumps(
            {
                "N": self.N,
                "clouds": [cloud.toJSON() for cloud in self.clouds],
                "labels": self.labels,
            }
        )

    def fromJSON(self, Jstring):
        J = json.loads(Jstring)
        self.N = J["N"]
        self.clouds = [FiducialCloud().fromJSON(js) for js in J["clouds"]]
        self.labels = J["labels"]
        return self

    def has_data(self):
        return all([fc.has_data() for fc in self.clouds])

    def data(self, idx=None, label=None):
        if not self.has_data():
            logger.warning("Data not loaded, cannot retrieve data")
            return

        if not (idx or label):
            raise ValueError("must provide either idx or label to retrieve data")
        elif label and not idx:
            try:
                idx = self.labels.index(label)
            except ValueError:
                raise ValueError(
                    "Could not find label {} in reg list: {}".format(label, self.labels)
                )
        elif label and not idx:
            print("Both label and index provided, using idx")
        return self.clouds[idx].data

    @property
    def count(self):
        return [c.count for c in self.clouds]

    @property
    def count_matching(self):
        return self.matching()[0].shape[1]

    @property
    def mincount(self):
        return [c.mincount for c in self.clouds]

    @mincount.setter
    def mincount(self, value):
        for c in self.clouds:
            c.mincount = value
        self._matching = self._get_matching()

    def _get_matching(self, inworld=False):
        """enforce matching points in cloudset"""
        if inworld:
            coords = [C.coords_inworld for C in self.clouds]
        else:
            coords = [C.coords for C in self.clouds]
        while True:
            # for n in range(1, self.N):
            #   self.clouds[0], self.clouds[n] = get_matching_points(
            #       self.clouds[0], self.clouds[n])
            # NOTE BY ATM
            # self.N == len(data): number of channels.
            # I think this goes through all combinations of pairs of
            # cloud points, 1 for each channel.
            #
            # For each pair, it finds sets of matching points
            # by comparing their XYZ coordinates and checking if
            # they fall within a maximum distance.
            # NOTE: is this in intrinsic, or world, space?
            # Maybe it's fine to do it intrinsic, because the answer
            # is the same either way.
            #
            # Then, such coordinates can be fed into a transformation,
            # such as scipy.ndimage.interpolation.map_coordinates,
            # to map the coordinates from a second set back to a first set.
            # This transformation is handled by imwarp(), and is invoked by
            # cloudset's tform(), like so:
            #
            # matching = self._get_matching(inworld=inworld)
            # moving   = matching[movIdx]
            # fixed    = matching[fixIdx]
            # tform    = funcDict[mode](moving, fixed)
            #
            # Somewhere (not clear?) the coordinates are converted
            # from "intrinsic" (pixel) to "world" (spatial), which I presume
            # is necessary to handle different pixel dimensions in
            # X/Y and in Z (making voxels which aren't "cubes").
            #
            # The different dx, dz values *only* matter when
            # fitting Gaussians, or when imwarping.
            #
            for m, n in itertools.combinations(range(self.N), 2):
                coords[m], coords[n] = get_matching_points(coords[m], coords[n])
            if len({c.shape for c in coords}) == 1:
                break
        if not all([len(c) for c in coords]):
            raise IndexError("At least one point cloud has no points")
        return coords

    def matching(self):
        if "_matching" in dir(self):
            return self._matching
        else:
            self._matching = self._get_matching()
            return self._matching

    def __getitem__(self, key):
        if isinstance(key, str) or (isinstance(key, int) and key > self.N):
            if self.labels is not None:
                if key in self.labels:
                    return self.clouds[self.labels.index(key)]
                else:
                    raise ValueError("Unrecognized label for CloudSet")
            else:
                raise ValueError(
                    "Cannot index CloudSet by string without " "provided labels"
                )
        elif isinstance(key, int) and key < self.N:
            return self.clouds[key]
        else:
            raise ValueError("Index must either be label or int < numClouds in Set")

    def get_all_tforms(
        self,
        refs=None,
        inworld=True,
        modes=(
            "translation",
            "rigid",
            "similarity",
            "affine",
            "2step",
            "cpd_2step",
            "cpd_similarity",
        ),
    ):
        """Generate an array of dicts for lots of possible tforms."""

        if refs is None:
            # default to all channels
            regto = self.labels
        else:
            # otherwise validate the provided list of channels
            try:
                iter(refs)
            except TypeError:
                refs = [refs]
            regto = []
            for ref in refs:
                if ref in self.labels:
                    regto.append(ref)
                else:
                    logger.warning(
                        "Reference {} not in lablels: {} ... skipping".format(
                            ref, self.labels
                        )
                    )
        if not len(regto):
            logger.error("No recognized values in refs list.  No tforms calculated")
            return None

        # validate modes, and assert iterable
        try:
            iter(modes)
        except TypeError:
            modes = [modes]
        modes = [m for m in modes if m in funcDict]

        import itertools

        pairings = itertools.permutations(self.labels, 2)
        pairings = [i for i in pairings if i[1] in regto]

        D = []
        for moving, fixed in pairings:
            for mode in modes:
                try:
                    D.append(
                        {
                            "mode": mode,
                            "reference": fixed,
                            "moving": moving,
                            "inworld": inworld,
                            "tform": self.tform(moving, fixed, mode, inworld=inworld),
                        }
                    )
                except Exception:
                    print("SKIPPING MODE: ", mode)
                    logger.error(
                        'Failed to calculate mode "{}" in get_all_tforms.  Skipping.'.format(
                            mode
                        )
                    )
        return D

    def get_all_tforms_dataframe(
        self,
        refs=None,
        inworld=True,
        modes=(
            "translation",
            "rigid",
            "similarity",
            "affine",
            "2step",
            "cpd_2step",
            "cpd_similarity",
        ),
    ):
        """
        Generate a dataframe for lots of possible tforms.
        Added by Andrew T. Martens.
        """

        if refs is None:
            # default to all channels
            regto = self.labels
        else:
            # otherwise validate the provided list of channels
            try:
                iter(refs)
            except TypeError:
                refs = [refs]
            regto = []
            for ref in refs:
                if ref in self.labels:
                    regto.append(ref)
                else:
                    logger.warning(
                        "Reference {} not in lablels: {} ... skipping".format(
                            ref, self.labels
                        )
                    )
        if not len(regto):
            logger.error("No recognized values in refs list.  No tforms calculated")
            return None

        # validate modes, and assert iterable
        try:
            iter(modes)
        except TypeError:
            modes = [modes]
        modes = [m for m in modes if m in funcDict]

        ### ATM EDIT
        ### From here, change how the data are stored & returned
        # NOTE for now, storing results in python lists,
        # and letting pandas convert to a dataframe.
        # Possibly could get a slight speed boost if
        # pre-allocating ndarrays.
        # However, issues with determining string sizes
        # & array sizes in advance. Could negate any speed gains.

        # Columns:
        # 1. mode (str)
        # 2. reference (str)
        # 3. moving (str)
        # 4. inworld (bool)
        # 5. tform (ndarray)
        #
        # Initialize numpy arrays with known size
        # size = len(modes) * len(regto) * len(regto) - 1

        modes_arr = []  # numpy.array(shape=size, dtype=str)
        refs_arr = []  # numpy.array(shape=size, dtype=str)
        moving_arr = []  # numpy.array(shape=size, dtype=str)
        inworld_arr = []  # numpy.array(shape=size, dtype=numpy.bool)
        tform_arr = []  # numpy.array(shape=(size, 4, 4), dtype=numpy.float64)

        pairings = itertools.permutations(self.labels, 2)
        pairings = [i for i in pairings if i[1] in regto]

        # If using ndarrays, then add entries by index
        # i = 0
        for moving, fixed in pairings:
            for mode in modes:
                try:
                    # NOTE another problem with ndarrays for str
                    # is to know what size to allocate in advance...
                    # modes_arr[i]  = mode
                    # refs_arr[i]   = fixed
                    # moving_arr[i] = moving
                    # inworld_arr   = inworld
                    # tform_arr[i]  = self.tform(moving, fixed, mode, inworld=inworld)

                    modes_arr.append(mode)
                    refs_arr.append(fixed)
                    moving_arr.append(moving)
                    inworld_arr.append(inworld)
                    tform_arr.append(self.tform(moving, fixed, mode, inworld=inworld))

                    # i += 1

                # NOTE that skipping will cause ndarrays to be longer than necessary...
                # For now, just use lists instead
                except Exception:
                    print("SKIPPING MODE: ", mode)
                    logger.error(
                        'Failed to calculate mode "{}" in get_all_tforms.  Skipping.'.format(
                            mode
                        )
                    )

        # Store the lists into a dict
        D = {
            "mode": modes_arr,
            "reference": refs_arr,
            "moving": moving_arr,
            "inworld": inworld_arr,
            "tform": tform_arr,
        }

        # Convert the dict of lists into a dataframe
        D = pandas.DataFrame(D)

        return D

    def write_all_tforms(self, outfile, **kwargs):
        """write all of the tforms for this cloudset to file"""

        class npEncoder(json.JSONEncoder):
            def fixedString(self, obj):
                numel = len(obj)
                form = "[" + ",".join(["{:14.10f}"] * numel) + "]"
                return form.format(*obj)

            def default(self, obj):
                if isinstance(obj, np.ndarray):
                    if all(isinstance(i, np.ndarray) for i in obj):
                        nestedList = obj.tolist()
                        result = [self.fixedString(l) for l in nestedList]
                        return result
                    else:
                        return obj.tolist()
                return json.JSONEncoder.default(self, obj)

        tforms = self.get_all_tforms(**kwargs)
        outdict = {"cloud": self.toJSON(), "tforms": tforms}
        outstring = json.dumps(outdict, cls=npEncoder, indent=2)
        outstring = outstring.replace('"[', " [").replace(']"', "]")
        with open(outfile, "w") as file:
            file.write(outstring)

    # Main Method
    def tform(
        self, movingLabel=None, fixedLabel=None, mode="2step", inworld=True, **kwargs
    ):
        """get tform matrix that maps moving point cloud to fixed point cloud

        Args:
            movingLabel (:obj:`str`): label/wave to register.  if none, will use the
                second cloud in the cloudset by default.
            fixedLabel (:obj:`str`):  reference label/wave.  if none, will use the
                first cloud in the cloudset by default.
            mode ({'translation', 'rigid', 'similarity', 'affine', '2step',
                'cpd_rigid', 'cpd_similarity', 'cpd_affine', 'cpd_2step'}):
                type of registration transformation to calculate.  CPD modes
                use coherent point drift estimation, other modes use least
                squares.
            inworld (:obj:`bool`): if True, will use :obj:`intrinsicToWorld` to
                convert pixel coordinates into world coordinates using the voxel
                size provided.  (needs work)

        Returns:
            :obj:`np.ndarray`: 4x4 transformation matrix


        """
        if self.labels is None:
            logging.warning("No label list provided... cannot get tform by label")
            return
        movingLabel = movingLabel if movingLabel is not None else self.labels[1]
        fixedLabel = fixedLabel if fixedLabel is not None else self.labels[0]

        try:
            movIdx = self.labels.index(movingLabel)
        except ValueError:
            raise ValueError(
                "Could not find label {} in reg list: {}".format(
                    movingLabel, self.labels
                )
            )
        try:
            fixIdx = self.labels.index(fixedLabel)
        except ValueError:
            raise ValueError(
                "Could not find label {} in reg list: {}".format(
                    fixedLabel, self.labels
                )
            )

        mode = mode.lower()
        if mode in funcDict:
            if mode.startswith("cpd"):
                if inworld:
                    moving = self.clouds[movIdx].coords_inworld.T
                    fixed = self.clouds[fixIdx].coords_inworld.T
                else:
                    moving = self.clouds[movIdx].coords.T
                    fixed = self.clouds[fixIdx].coords.T
                if "2step" in mode:
                    tform = funcDict[mode](moving, fixed)
                else:
                    reg = funcDict[mode](moving, fixed)
                    tform = reg.register(None)[4]
            else:
                matching = self._get_matching(inworld=inworld)
                moving = matching[movIdx]
                fixed = matching[fixIdx]
                tform = funcDict[mode](moving, fixed)

                if inworld:
                    # FIXME: Why the hell is this necessary??
                    # pointed out by Rainer... but after looking through the code for a while,
                    # and carefully checking all the intrinsic-world conversions,
                    # I still can't figure out why the tforms are all off by half of a pixel
                    tform[0][3] -= self.dx / 2
                    tform[1][3] -= self.dx / 2
                    tform[2][3] -= self.dz / 2
            logger.info(
                "Measured {} Tform Matrix {}inWorld:\n".format(
                    mode, "" if inworld else "not "
                )
                + str(tform)
            )
            return tform
        else:
            raise ValueError("Unrecognized transformation mode: {}".format(mode))

    def show(self, matching=False, withimage=True, filtered=True):
        """show points in clouds overlaying image, if matching is true, only
        show matching points from all sets"""
        plt.figure(figsize=(12, 12))
        if withimage:
            if not self.has_data():
                raise AttributeError(
                    "No data present in cloudset.  Cannot show"
                    " image.  Try calling reload_data() from parent object"
                )
            if filtered and self.clouds[0].filtered is not None:
                im = self.clouds[0].filtered.max(0)
            else:
                im = self.clouds[0].max(0)
            plt.imshow(im, cmap="gray", vmax=im.max() * 0.7)

        colors = ["red", "purple", "magenta", "blue", "green"]
        for i in reversed(range(self.N)):
            if matching:
                X = self.matching()[i][0]
                Y = self.matching()[i][1]
            else:
                X = self.clouds[i].coords[0]
                Y = self.clouds[i].coords[1]
            plt.scatter(X, Y, c=colors[i], s=5)
        plt.show()

    def show_matching(self, **kwargs):
        self.show(matching=True, **kwargs)

    def show_tformed(self, movingLabel=None, fixedLabel=None, matching=False, **kwargs):
        T = self.tform(movingLabel, fixedLabel, **kwargs)
        if matching:
            movingpoints = self.matching()[self.labels.index(movingLabel)]
            fixedpoints = self.matching()[self.labels.index(fixedLabel)]

            matching = self._get_matching(inworld=kwargs.get("inworld", False))
            movingpoints = matching[self.labels.index(movingLabel)]
            fixedpoints = matching[self.labels.index(fixedLabel)]

        else:
            movingpoints = self[movingLabel].coords
            fixedpoints = self[fixedLabel].coords
        shiftedpoints = affineXF(movingpoints, T)
        fp = plt.scatter(fixedpoints[0], fixedpoints[1], c="b", s=7)
        mp = plt.scatter(movingpoints[0], movingpoints[1], c="m", marker="x", s=5)
        sp = plt.scatter(shiftedpoints[0], shiftedpoints[1], c="r", s=7)
        plt.legend((fp, mp, sp), ("Fixed", "Moving", "Registered"))
        plt.show()

    def show_tformed_matching(self, movingLabel=None, fixedLabel=None, **kwargs):
        self.show_tformed(movingLabel, fixedLabel, True, **kwargs)

    def show_tformed_image(self, movingLabel=None, fixedLabel=None, **kwargs):
        try:
            from llspy.libcudawrapper import affineGPU
        except ImportError:
            print("Could not import affineGPU, can't show tformed image")
            return
        T = self.tform(movingLabel, fixedLabel, **kwargs)
        movingImg = self.data(label=movingLabel)
        fixedImg = self.data(label=fixedLabel)
        movingReg = affineGPU(movingImg, T, [self.dz, self.dx, self.dx])
        imshowpair(fixedImg, movingReg, **kwargs)


def imoverlay(im1, im2, method=None, mip=False):
    im1 = im1.astype(np.float) if not mip else im1.astype(np.float).max(0)
    im2 = im2.astype(np.float) if not mip else im2.astype(np.float).max(0)
    im1 -= im1.min()
    im1 /= im1.max()
    im2 -= im2.min()
    im2 /= im2.max()

    ndim = im1.ndim
    if method == "diff":
        im3 = im1 - im2
        im3 -= im3.min()
        im3 /= im3.max()
        return im3
    else:
        return np.stack((im1, im2, im1), ndim)


def imshowpair(im1, im2, method=None, mip=False, **kwargs):
    # normalize
    if not im1.shape == im2.shape:
        raise ValueError("images must be same shape")

    if not mip:
        try:
            from tifffile import imshow
        except ImportError:
            imshow = plt.imshow
            mip = True
    else:
        imshow = plt.imshow
        if im1.ndim < 3:
            mip = False

    if method == "diff":
        imshow(imoverlay(im1, im2, "diff"), cmap="gray", vmin=0.2, vmax=0.8)
    elif method == "3D":
        im3 = imoverlay(im1, im2)
        fig, subpl, ax = imshow(im3, subplot=221)
        imshow(np.rot90(im3.max(1)), figure=fig, subplot=222)
        imshow(im3.max(2), figure=fig, subplot=223)
    else:  # falsecolor
        imshow(imoverlay(im1, im2))

    plt.show()


class RegFile(object):
    def __init__(self, path):
        self.path = path
        if not osp.isfile(path):
            raise FileNotFoundError("Could not find registration file: {}".format(path))
        self.parsefile()

    def parsefile(self):
        try:
            with open(self.path) as json_data:
                regdict = json.load(json_data)
        except json.decoder.JSONDecodeError:
            raise

        if "tforms" not in regdict:
            self.tforms = []
            return
        else:
            self.tforms = regdict["tforms"]

        # these parameters are written to the regfile by llspy
        try:
            from datetime import datetime

            self.date = regdict.get("date", None)
            if isinstance(self.date, list):  # backwards compatibility
                self.date = self.date[0]
            if self.date:
                self.date = datetime.strptime(self.date, "%Y/%m/%d-%H:%M")
        except Exception as e:
            logger.error("Could not parse registration file date: {}".format(e))

        self.dx = regdict.get("dx", None)
        self.dz = regdict.get("dz", None)
        self.z_motion = regdict.get("z_motion", None)
        self.regdir = regdict.get("path", None)

        self.tform_dict = {}
        self.refwaves = []
        self.movwaves = []
        self.modes = []
        for tform in self.tforms:
            ref = str(tform["reference"])
            mov = str(tform["moving"])
            mode = tform["mode"].lower()
            self.refwaves.append(ref)
            self.movwaves.append(mov)
            self.modes.append(mode)
            if ref not in self.tform_dict:
                self.tform_dict[ref] = {}
            if mov not in self.tform_dict[ref]:
                self.tform_dict[ref][mov] = {}
            self.tform_dict[ref][mov][mode] = tform["tform"]
        self.refwaves = sorted(list(set(self.refwaves)))
        self.movwaves = sorted(list(set(self.movwaves)))
        self.modes = sorted(list(set(self.modes)))
        self.waves = self.refwaves  # to make it easier to substitute for RegDir

    @property
    def n_tforms(self):
        return len(self.tforms)

    @property
    def isValid(self):
        return bool(len(self.tforms))

    def get_tform(self, moving, ref, mode):
        ref = str(ref)
        moving = str(moving)
        mode = str(mode).lower()
        if ref not in self.tform_dict:
            raise RegistrationError(
                "Reference wave {} not in registration file".format(ref)
            )
        if moving not in self.tform_dict[ref]:
            raise RegistrationError(
                "No transform to map moving wave {} onto refrence wave {}".format(
                    moving, ref
                )
            )
        if mode not in self.tform_dict[ref][moving]:
            raise RegistrationError(
                "Transform mode {} not found for refwave: {}, movingwave: {}".format(
                    mode, ref, moving
                )
            )

        return self.tform_dict[ref][moving][mode]


funcDict = {
    "translate": infer_translation,
    "translation": infer_translation,
    "rigid": infer_rigid,
    "similarity": infer_similarity,
    "affine": infer_affine,
    "2step": infer_2step,
    "cpd_rigid": CPDrigid,
    "cpd_similarity": CPDsimilarity,
    "cpd_affine": CPDaffine,
    "cpd_2step": cpd_2step,
}
