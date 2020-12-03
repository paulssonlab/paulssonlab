#!/usr/bin/env python

from magicgui import magicgui
import napari
from napari.layers import Image, Labels

##

import tables
import numpy

from algorithms import niblack_phase_segmentation
from trench_detect import find_trenches
from params import read_params_file
from metadata import get_metadata

from scipy.signal import find_peaks

###


def detect_rows(img, axis, cutoff, plateau_size):
    """
    TODO I think it might be worth adding a manual padding ability,
    in case it's just slightly off. Or a a way to constrain or define
    the final dimension (length of trenches).
    """
    # The variance is expected to be greater where there are features (numerals or trenches)
    var = numpy.var(img, axis=axis)

    # The mean is useful as a normalization
    mean = numpy.mean(img, axis=axis)

    # Normalize the variance by the mean
    var_mean = var / mean

    # Define an arbitrary threshold, below which we ignore,
    # and above which we deem peak-worthy
    gt = var_mean > cutoff

    # By defining the plateau size as such, we can eliminate the peaks
    # which arise from the emblazoned numerals.
    peaks, peaks_properties = find_peaks(gt, plateau_size=plateau_size)

    # FIXME img.shape[0] or [1] ? (width, or height?)
    rows = numpy.array(
        [
            [0, i[0], img.shape[1], i[1]]
            for i in zip(
                peaks_properties["left_edges"], peaks_properties["right_edges"]
            )
        ]
    )

    return rows


def regions_arr_to_canvas(canvas, array):
    """
    Input a numpy array of regions, and a canvas (which may or may not be empty).
    Return a canvas containing the labeled regions.
    """
    for i, r in enumerate(array):
        canvas[r[0] : r[2], r[1] : r[3]] = i + 1

    return canvas


###

### Begin main

# Load a sample microscopy image from HDF5
in_file_path = (
    "/home/andrew/Documents/Harvard Work/Microscopy Data/2020-11-25/HDF5/data.h5"
)
h5file = tables.open_file(in_file_path, "r")
FILE = "File_run"
FOV = "FOV_0"
FRAME = "Frame_0"
Z = "Z_0"

metadata_node = h5file.get_node("/Metadata/{}".format(FILE))()
CHANNELS = get_metadata(metadata_node)["channels"]
px_mu = get_metadata(metadata_node)["pixel_microns"]
# NOTE could add aditional information, such as a mapping between channels and colormaps in the viewer.


# Load the params from a file
params_trenches_and_rows = read_params_file("trench_params.yaml")
trench_row_params = params_trenches_and_rows["trench_rows"]
trench_params = params_trenches_and_rows["trenches"]

params_segmentation = read_params_file("seg_params.yaml")


###

# NOTE: Define the GUIs inside the main function so that we can load the params from a file
# This feels really clunky, defining a function within a function.
# For one, it prohibits re-use outside of this function.
# Is there no other way?

### Make the Trench Row GUI


@magicgui(
    auto_call=False,
    layout="vertical",
    call_button="Calculate",
    plateau_size={"minimum": 0, "maximum": 1000},
    cutoff={"minimum": 0, "maximum": 100},
)
def run_trench_row_detection_napari(
    Layer_Phase: Image,
    plateau_size: float = trench_row_params["plateau_size"],
    cutoff: int = trench_row_params["cutoff"],
) -> Labels:

    kwargs = {
        "axis": 0,
        "plateau_size": int(plateau_size / px_mu),
        "cutoff": cutoff,
    }

    img = Layer_Phase.data
    img = numpy.rot90(img, 1)

    # Get a numpy array of regions
    ranges = detect_rows(img, **kwargs)

    # Make a labeled image
    # FIXME or height, width?
    canvas = numpy.zeros(shape=(img.shape[0], img.shape[1]), dtype=numpy.uint16)
    canvas = regions_arr_to_canvas(canvas, ranges)
    canvas = numpy.rot90(canvas, -1)

    # TODO how to pass in a name for the resulting labels layer?
    return canvas


### Make the Trench GUI


@magicgui(
    auto_call=False,
    layout="vertical",
    call_button="Calculate",
    unsharp_mask_radius={"minimum": 0, "maximum": 10},
    unsharp_mask_amount={"minimum": 0, "maximum": 100},
    distance={"minimum": 0, "maximum": 1000},
    plateau_size={"minimum": 0, "maximum": 1000},
    cutoff={"minimum": 0.0, "maximum": 1.0},
)
def run_trench_detection_napari(
    Layer_Phase: Image,
    # Modifiable values for trench row detection:
    unsharp_mask_radius: int = trench_params["unsharp_mask_radius"],
    unsharp_mask_amount: float = trench_params["unsharp_mask_amount"],
    distance: float = trench_params["distance"],
    plateau_size: float = trench_params["plateau_size"],
    cutoff: float = trench_params["cutoff"],
    # TODO pass in microns, convert to pixels
    crop_top: int = trench_params["crop_top"],
    crop_bottom: int = trench_params["crop_bottom"],
) -> Labels:
    # Temporarily assign values to any variables which aren't connected to magicgui
    kwargs = {
        "unsharp_mask_radius": unsharp_mask_radius,
        "unsharp_mask_amount": unsharp_mask_amount,
        "axis": 1,
        "distance": int(distance / px_mu),
        "plateau_size": int(plateau_size / px_mu),
        "cutoff": cutoff,
        "crop_top": crop_top,
        "crop_bottom": crop_bottom,
    }

    # FIXME is there a way to pass in the rows which were detected using the previous GUI widgets?
    rows = numpy.array([[0, 383, 2048, 553], [0, 1181, 2048, 1335]])

    img = Layer_Phase.data

    # FIXME or height, width?
    canvas = numpy.zeros(shape=(img.shape[0], img.shape[1]), dtype=numpy.uint16)

    for row in rows:
        # Rotate
        # NOTE Apparently all Napari layers are loaded as F32. Doing any computation
        # might vary if expecting uint16!
        img = img.astype(numpy.uint16)
        img = numpy.rot90(img, 1)
        img_slice = img[row[0] : row[2], row[1] : row[3]]

        regions = find_trenches(img_slice, **kwargs)

        # Fix the coordinates to match the whole image
        # NOTE it's also possible to create a smaller canvas,
        # and to apply the "translate" parameter to it in Napari.
        # Could be better? Would that make the coordinates confusing?
        # Use less memory?
        for r in regions:
            r[0] += row[0]
            r[1] += row[1]
            r[2] += row[0]
            r[3] += row[1]

        # Make a canvas
        canvas = regions_arr_to_canvas(canvas, regions)

    canvas = numpy.rot90(canvas, -1)

    return canvas


@magicgui(
    auto_call=False,
    layout="vertical",
    call_button="Calculate",
    niblack_k={"minimum": -1.0, "maximum": 1.0},
    # Horizontal window size
    # Must be an odd number! TODO How to force or enforce this constraint?
    niblack_w_h={"minimum": 1, "maximum": 21},
    # Vertical window size
    # Must be an odd number! TODO How to force or enforce this constraint?
    niblack_w_v={"minimum": 1, "maximum": 21},
    otsu_multiplier={"maximum": 10.0},
    fluor_background={"minimum": 0.0, "maximum": 65536.0},
    garbage_otsu_value={"minimum": 0, "maximum": 65536},
    scaling_factor={"minimum": 0.0, "maximum": 10.0},
    fluor_sigma={"minimum": 0.0, "maximum": 10.0},
    phase_background={"minimum": 0.0, "maximum": 10.0},
    phase_threshold_min={"minimum": 0.0, "maximum": 65536.0},
    phase_threshold_max={"minimum": 0.0, "maximum": 65536.0},
    phase_sigma={"minimum": 0.0, "maximum": 10.0},
    min_size={"minimum": 0, "maximum": 1000},
)
def run_niblack_phase_segmentation_napari(
    Layer_Fluor: Image,
    Layer_Phase: Image,
    niblack_k: float = 0.0,
    niblack_w_h: int = 7,
    niblack_w_v: int = 7,
    otsu_multiplier: float = 1.0,
    fluor_background: float = 0.0,
    garbage_otsu_value: int = 200,
    scaling_factor: float = 1.0,
    fluor_sigma: float = 1.0,
    phase_background: float = 0.0,
    phase_threshold_min: float = 0.0,
    phase_threshold_max: float = 65536.0,
    phase_sigma: float = 1.9,
    min_size: float = 9.75,
) -> Labels:
    # Take the 2 layers & convert each into an appropriate stack-like object
    # Initalize the empty ndarray
    stack = numpy.empty(
        (
            Layer_Fluor.data.shape[0],
            Layer_Fluor.data.shape[1],
            1,
            2,
        ),  # 1 region, 2 channels
        dtype=numpy.float32,
        order="F",
    )

    # TODO change these so that, instead of loading the entire images,
    # it uses detected trench rows or trenches.

    stack[..., 0, 0] = Layer_Fluor.data
    stack[..., 0, 1] = Layer_Phase.data

    stack_fl = stack[..., 0]
    stack_ph = stack[..., 1]

    # Temporarily assign values to any variables which aren't connected to magicgui
    kwargs = {
        "niblack_k": niblack_k,
        "niblack_w": [niblack_w_h, niblack_w_v],
        "otsu_multiplier": otsu_multiplier,
        "fluor_background": fluor_background,
        "garbage_otsu_value": garbage_otsu_value,
        "scaling_factor": scaling_factor,
        "fluor_sigma": fluor_sigma,
        "phase_background": phase_background,
        "phase_threshold_min": phase_threshold_min,
        "phase_threshold_max": phase_threshold_max,
        "phase_sigma": phase_sigma,
        "min_size": int(min_size / px_mu),
    }

    # Call the function & get a labeled image
    labeled = niblack_phase_segmentation(stack_fl, stack_ph, **kwargs)
    return labeled[..., 0]


# Main GUI init:
with napari.gui_qt():
    # create a viewer and add some images
    viewer = napari.Viewer()

    for c in CHANNELS:
        c = c.decode("utf-8")
        img = h5file.get_node(
            "/Images/{}/{}/{}/{}/{}".format(FILE, FOV, FRAME, Z, c)
        ).read()
        img = numpy.rot90(img, -1)
        # NOTE it's possible to add rotate=-90 as one of the parameters,
        # however this causes the coordinates to become negative :(
        # Unclear how to do a rotation and keep normal coords.
        viewer.add_image(img.astype("float"), name=c)

    # Trench Rows

    # instantiate the widget
    gui_detect_rows = run_trench_row_detection_napari.Gui()
    gui_detect_rows.result_name = "Detected Rows"

    # add the gui to the viewer as a dock widget
    # TODO is there also a way to manually add a widget (the first one) which is just a Label?
    viewer.window.add_dock_widget(gui_detect_rows, name="Trench Rows", area="left")

    # if a layer gets added or removed, refresh the dropdown choices
    # NOTE why Layer_Phase? Does that set the default?
    viewer.layers.events.changed.connect(
        lambda x: gui_detect_rows.refresh_choices("Layer_Phase")
    )

    # TODO Assume that at this point we will have access to the rows which were returned by the previous widget.

    # Trenches

    # instantiate the widget
    gui_detect_trenches = run_trench_detection_napari.Gui()
    gui_detect_trenches.result_name = "Detected Trenches"

    # add the gui to the viewer as a dock widget
    viewer.window.add_dock_widget(gui_detect_trenches, name="Trenches", area="left")

    # if a layer gets added or removed, refresh the dropdown choices
    viewer.layers.events.changed.connect(
        lambda x: gui_detect_trenches.refresh_choices("Layer_Phase")
    )

    # Segmentation

    # instantiate the widget
    gui_detect_trenches = run_niblack_phase_segmentation_napari.Gui()
    gui_detect_trenches.result_name = "Segmented Cells"

    # add the gui to the viewer as a dock widget
    viewer.window.add_dock_widget(gui_detect_trenches, name="Segmentation", area="left")

    # if a layer gets added or removed, refresh the dropdown choices
    viewer.layers.events.changed.connect(
        lambda x: gui_detect_trenches.refresh_choices("Layer_Phase")
    )

h5file.close()
