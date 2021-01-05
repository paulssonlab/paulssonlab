#!/usr/bin/env python

import click

##

from magicgui import magicgui
import napari
from napari.layers import Image, Labels

##

import tables
import numpy

from algorithms import niblack_phase_segmentation
from params import read_params_file
from metadata import get_metadata

from new_trench_algos import find_trenches, generate_boxes
from napari_browse_hdf5 import (
    add_image_layers,
    get_largest_extents_hdf5,
    metadata_attributes_equal,
    metadata_array_equal,
)

from scipy.signal import find_peaks

###

# def detect_rows(img, axis, cutoff, plateau_size):
# """
# TODO I think it might be worth adding a manual padding ability,
# in case it's just slightly off. Or a a way to constrain or define
# the final dimension (length of trenches).
# """
## The variance is expected to be greater where there are features (numerals or trenches)
# var = numpy.var(img, axis=axis)

## The mean is useful as a normalization
# mean = numpy.mean(img, axis=axis)

## Normalize the variance by the mean
# var_mean = var/mean

## Define an arbitrary threshold, below which we ignore,
## and above which we deem peak-worthy
# gt = var_mean > cutoff

## By defining the plateau size as such, we can eliminate the peaks
## which arise from the emblazoned numerals.
# peaks, peaks_properties = find_peaks(gt, plateau_size=plateau_size)

## FIXME img.shape[0] or [1] ? (width, or height?)
# rows = numpy.array([ [0, i[0], img.shape[1], i[1]] for i in zip (peaks_properties["left_edges"], peaks_properties["right_edges"]) ])

# return rows


def regions_arr_to_canvas(canvas, array):
    """
    Input a numpy array of regions, and a canvas (which may or may not be empty).
    Return a canvas containing the labeled regions.
    """
    for i, r in enumerate(array):
        canvas[r[0] : r[2], r[1] : r[3]] = i + 1

    return canvas


def regions_arr_to_canvas_F_order(canvas, regions, labels):
    """
    Input a numpy array of regions, and a canvas (which may or may not be empty).
    Return a canvas containing the labeled regions.
    For F-ordered arrays.
    """
    for i, r in enumerate(regions):
        canvas[r[0] : r[2], r[1] : r[3]] = labels[..., i]

    # for i in range(array.shape[2]):
    # r = array[..., i]
    # canvas[r[0] : r[2], r[1] : r[3]] = i + 1

    return canvas


###

###


def main_function(
    in_file,
    FILE,
    FOV,
    FRAME,
    Z_LEVEL,
    napari_settings_file,
    params_file_trenches,
    params_file_segmentation,
):
    """
    This function should be called by the Click command-line interface. Pass in:
    - path to hdf5 file containing images
    - desired file, fov, frame, z_level to open for parameter testing
    - yaml file for how to load images
    - yaml file for trench detection params
    - yaml file for segmentation after trench detection

    NOTE right now we don't use the napari settings yet??
    """

    # Will this allow for the trenches calculated by trench detection
    # to automatically become available to subsequent calculations?
    # A bit of a kludge, but good enough for now.
    global trenches
    trenches = None

    # Load a sample microscopy image from HDF5
    h5file = tables.open_file(in_file, "r")

    metadata_node = h5file.get_node("/Metadata/{}".format(FILE))()
    CHANNELS = get_metadata(metadata_node)["channels"]
    px_mu = get_metadata(metadata_node)["pixel_microns"]

    # TODO napari viewer settings
    # NOTE could add aditional information, such as a mapping between channels and colormaps in the viewer.
    layer_params = read_params_file(napari_settings_file)

    # Load the params from a file
    params_trenches_and_rows = read_params_file(params_file_trenches)
    # trench_row_params = params_trenches_and_rows["trench_rows"]
    trench_params = params_trenches_and_rows["trenches"]

    params_segmentation = read_params_file(params_file_segmentation)
    params_segmentation = params_segmentation["mKate2"]["parameters"]

    ###

    # NOTE: Define the GUIs inside the main function so that we can load the params from a file
    # This feels really clunky, defining a function within a function.
    # For one, it prohibits re-use outside of this function.
    # Is there no other way?

    #### Make the Trench Row GUI

    # @magicgui(
    # auto_call    = False,
    # layout       = "vertical",
    # call_button  = "Calculate",
    # plateau_size = {"minimum": 0,   "maximum": 1000},
    # cutoff       = {"minimum": 0, "maximum": 100},
    # )
    # def run_trench_row_detection_napari(
    # Layer_Phase: Image,
    # plateau_size : float = trench_row_params["plateau_size"],
    # cutoff       : int = trench_row_params["cutoff"],
    # ) -> Labels:

    # kwargs = {
    # "axis"         : 0,
    # "plateau_size" : int(plateau_size / px_mu),
    # "cutoff"       : cutoff,
    # }

    # img = Layer_Phase.data
    # img = numpy.rot90(img, 1)

    ## Get a numpy array of regions
    # ranges = detect_rows(img, **kwargs)

    ## Make a labeled image
    ## FIXME or height, width?
    # canvas = numpy.zeros(shape=(img.shape[0], img.shape[1]), dtype=numpy.uint16)
    # canvas = regions_arr_to_canvas(canvas, ranges)
    # canvas = numpy.rot90(canvas, -1)

    ## TODO how to pass in a name for the resulting labels layer?
    # return canvas

    ### Make the Trench GUI

    @magicgui(
        auto_call=False,
        layout="vertical",
        call_button="Calculate",
        tr_height={"minimum": 0, "maximum": 100},
        tr_width={"minimum": 0, "maximum": 100},
        tr_spacing={"minimum": 0, "maximum": 1000, "decimals": 3},
        pad_left={"minimum": 0, "maximum": 10},
        pad_right={"minimum": 0, "maximum": 10},
        pad_top={"minimum": 0, "maximum": 20},
        pad_bottom={"minimum": 0, "maximum": 20},
    )
    def run_trench_detection_napari(
        Layer_Phase: Image,
        tr_height: float = trench_params["tr_height"],
        # tr_width:   float = trench_params["tr_width"],
        tr_width: int = trench_params["tr_width"],
        tr_spacing: float = trench_params["tr_spacing"],
        pad_left: int = trench_params["pad_left"],
        pad_right: int = trench_params["pad_right"],
        pad_top: int = trench_params["pad_top"],
        pad_bottom: int = trench_params["pad_bottom"],
    ) -> Labels:
        kwargs = {
            "tr_height": int(tr_height / px_mu),
            # "tr_width"   : int(tr_width / px_mu),
            "tr_width": tr_width,
            # "tr_spacing" : tr_spacing / px_mu,
            "tr_spacing": tr_spacing,
            "pad_left": pad_left,
            "pad_top": pad_top,
            "pad_right": pad_right,
            "pad_bottom": pad_bottom,
        }

        img = Layer_Phase.data

        # find_trenches() returns a multi-dimensional numpy list (and not a numpy array,
        # because each row might have a different number of trenches,
        # after taking into account out-of-bounds issues).
        boxes = generate_boxes(
            kwargs["tr_width"], kwargs["tr_height"], kwargs["tr_spacing"], img.shape[0]
        )

        global trenches
        (trenches, y_offsets) = find_trenches(img.T, boxes, **kwargs)

        # FIXME or height, width?
        canvas = numpy.zeros(shape=(img.shape[0], img.shape[1]), dtype=numpy.uint16)

        # For now, share the same canvas
        # It would be nice if each row of trenches had its own canvas!
        for r in trenches:
            canvas = regions_arr_to_canvas(canvas, r)

        canvas = canvas.T

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
        niblack_w_v={"minimum": 1, "maximum": 41},
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
        # Layer_Trenches: Labels,
        niblack_k: float = params_segmentation["niblack_k"],
        # NOTE: assumes that the params are passed in as a list of 2 values
        niblack_w_h: int = params_segmentation["niblack_w"][0],
        niblack_w_v: int = params_segmentation["niblack_w"][1],
        otsu_multiplier: float = params_segmentation["otsu_multiplier"],
        fluor_background: float = params_segmentation["fluor_background"],
        garbage_otsu_value: int = params_segmentation["garbage_otsu_value"],
        scaling_factor: float = params_segmentation["scaling_factor"],
        fluor_sigma: float = params_segmentation["fluor_sigma"],
        phase_background: float = params_segmentation["phase_background"],
        phase_threshold_min: float = params_segmentation["phase_threshold_min"],
        phase_threshold_max: float = params_segmentation["phase_threshold_max"],
        phase_sigma: float = params_segmentation["phase_sigma"],
        min_size: float = params_segmentation["min_size"],
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
        # Must state that var. is global so that it can use the info from trench detection routine!
        global trenches
        # TODO change these so that, instead of loading the entire images,
        # it uses detected trench rows or trenches.

        # FIXME or height, width?
        canvas = numpy.zeros(
            shape=(Layer_Fluor.shape[0], Layer_Fluor.shape[1]), dtype=numpy.uint16
        )

        # First row of trenches, previously detected using another GUI
        # NOTE Here, regions is a *global* variable, previously
        # allocated by pressing calculate button when detecting trenches.
        for regions in trenches:
            ch_to_index = {}
            x_dimension = regions[0, 2] - regions[0, 0]
            y_dimension = regions[0, 3] - regions[0, 1]
            # Initialize the array as empty, with the 'F' ordering, and then to load it x/y/r/c
            stack = numpy.empty(
                (x_dimension, y_dimension, regions.shape[0], 2),
                dtype=numpy.float32,
                order="F",
            )

            # min_row, min_col, max_row, max_col
            # Load fluor data
            this_img = Layer_Fluor.data.T
            for j, r in enumerate(regions):
                # NOTE Does copy flag actually make a difference?
                stack[..., j, 0] = this_img[r[0] : r[2], r[1] : r[3]].astype(
                    numpy.float64, casting="safe", copy=False
                )

            # Load phase data
            this_img = Layer_Phase.data.T
            for j, r in enumerate(regions):
                # NOTE Does copy flag actually make a difference?
                stack[..., j, 1] = this_img[r[0] : r[2], r[1] : r[3]].astype(
                    numpy.float64, casting="safe", copy=False
                )

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

            stack_fl = stack[..., 0]
            stack_ph = stack[..., 1]

            # Call the function & get a labeled image
            labeled = niblack_phase_segmentation(stack_fl, stack_ph, **kwargs)

            canvas = regions_arr_to_canvas_F_order(canvas, regions, labeled)

        # Done writing to canvas
        canvas = canvas.T
        return canvas

    # Main GUI init:
    with napari.gui_qt():
        # create a viewer and add some images
        viewer = napari.Viewer()

        ## Get the largest extents across all the nd2 files within the h5 file,
        ## so that the dask array is made with the proper min, max dimensions for each dimension.
        # attributes = ["fields_of_view", "frames", "z_levels"]
        # extents = {}
        # for a in attributes:
        # (smallest, largest) = get_largest_extents_hdf5(h5file, a)
        # if smallest == largest:
        # extents[a] = [smallest]
        # else:
        # extents[a] = [i for i in range(smallest, largest + 1)]

        ## Get the height & width, while checking that they are identical across all nd2 files
        # height = metadata_attributes_equal(h5file, "height")
        # width = metadata_attributes_equal(h5file, "width")

        ## Define the channels
        # channels = metadata_array_equal(h5file, "channels")

        # corrections_file = None

        ## Image layers
        # add_image_layers(
        # h5file,
        # height,
        # width,
        # extents["fields_of_view"],
        # extents["frames"],
        # extents["z_levels"],
        # viewer,
        # channels,
        # layer_params,
        # corrections_file,
        # )

        for c in CHANNELS:
            c = c.decode("utf-8")
            img = h5file.get_node(
                "/Images/{}/FOV_{}/Frame_{}/Z_{}/{}".format(
                    FILE, FOV, FRAME, Z_LEVEL, c
                )
            ).read()
            img = img.T
            # NOTE it's possible to add rotate=-90 as one of the parameters,
            # however this causes the coordinates to become negative :(
            # Unclear how to do a rotation and keep normal coords.
            viewer.add_image(img.astype("float"), **layer_params[c])

        ## Trench Rows

        ## instantiate the widget
        # gui_detect_rows = run_trench_row_detection_napari.Gui()
        # gui_detect_rows.result_name = "Detected Rows"

        ## add the gui to the viewer as a dock widget
        ## TODO is there also a way to manually add a widget (the first one) which is just a Label?
        # viewer.window.add_dock_widget(gui_detect_rows, name="Trench Rows", area="left")

        ## if a layer gets added or removed, refresh the dropdown choices
        ## NOTE why Layer_Phase? Does that set the default?
        # viewer.layers.events.changed.connect(lambda x: gui_detect_rows.refresh_choices("Layer_Phase"))

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
        viewer.window.add_dock_widget(
            gui_detect_trenches, name="Segmentation", area="left"
        )

        # if a layer gets added or removed, refresh the dropdown choices
        viewer.layers.events.changed.connect(
            lambda x: gui_detect_trenches.refresh_choices("Layer_Phase")
        )

    h5file.close()


### Begin main


@click.group()
def cli():
    """Invoke a Command to perform the desired operation:"""
    pass


@cli.command(no_args_is_help=True)
@click.option(
    "-i",
    "--in-file",
    "in_file",
    required=True,
    default="ND2/data.h5",
    type=str,
    help="Input HDF5 file.",
    show_default=True,
)
@click.option(
    "--napari-settings-file",
    "napari_settings_file",
    required=False,  # TODO implement this & set to True
    default="napari_settings.yaml",
    type=str,
    help="Napari settings file (YAML).",
    show_default=True,
)
@click.option(
    "--params-file-trenches",
    "params_file_trenches",
    required=True,
    default="trench_params.yaml",
    type=str,
    help="Regions detection parameters file (YAML).",
    show_default=True,
)
@click.option(
    "--params-file-segmentation",
    "params_file_segmentation",
    required=True,
    default="seg_params.yaml",
    type=str,
    help="Segmentation parameters file (YAML).",
    show_default=True,
)
@click.option(
    "--file",
    "file",
    required=True,
    default=None,
    type=str,
    help="File_ + name of original ND2 file",
)
@click.option(
    "--fov",
    "fov",
    required=True,
    default=0,
    type=int,
    help="FOV num",
    show_default=True,
)
@click.option(
    "--frame",
    "frame",
    required=True,
    default=0,
    type=int,
    help="Frame num",
    show_default=True,
)
@click.option(
    "--z_level",
    "z_level",
    required=True,
    default=0,
    type=int,
    help="z level",
    show_default=True,
)
def browse_parameters(
    in_file,
    file,
    fov,
    frame,
    z_level,
    napari_settings_file,
    params_file_trenches,
    params_file_segmentation,
):
    """Use Napari pick trench detection & cell segmentation parameters."""
    main_function(
        in_file,
        file,
        fov,
        frame,
        z_level,
        napari_settings_file,
        params_file_trenches,
        params_file_segmentation,
    )


###

if __name__ == "__main__":
    cli()

    # rows = find_trench_rows(numpy.rot90(img, 1), tr_width, int(tr_height / px_mu), tr_spacing)

    ## FIXME or height, width?
    # canvas = numpy.zeros(shape=(img.shape[0], img.shape[1]), dtype=numpy.uint16)

    # for row in rows:
    ## Rotate
    ## NOTE Apparently all Napari layers are loaded as F32. Doing any computation
    ## might vary if expecting uint16!
    # img = img.astype(numpy.uint16)
    # img = numpy.rot90(img, 1)
    # img_slice = img[row[0] : row[2], row[1] : row[3]]

    # trenches = find_trenches(img_slice)

    ## Fix the coordinates to match the whole image
    ## NOTE it's also possible to create a smaller canvas,
    ## and to apply the "translate" parameter to it in Napari.
    ## Could be better? Would that make the coordinates confusing?
    ## Use less memory?
    # for tr in trenches:
    # tr[0] += row[0]
    # tr[1] += row[1]
    # tr[2] += row[0]
    # tr[3] += row[1]

    ## Make a canvas
    # canvas = regions_arr_to_canvas(canvas, trenches)
