#!/usr/bin/env python

import click

# from dask import delayed
# import dask.array

from magicgui import magicgui
import napari
from napari.layers import Image, Labels

import tables
import numpy

from algo_sharp import fluor_sharpen_segmentation
from params import read_params_file
from metadata import get_metadata

from new_trench_algos import find_trenches, generate_boxes
from napari_browse_hdf5 import (
    add_image_layers,
    get_largest_extents_hdf5,
    metadata_attributes_equal,
    metadata_array_equal,
)

###


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

    return canvas


###


def edit_params(params):
    """
    Go through "parameters" and just write over its value.
    We can't do px_mu conversion, because the widgets also have to.
    """
    p = params["parameters"]["trenches"]
    for k, parameter in p.items():
        p[k] = p[k]["value"]

    return p


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
    trench_params = read_params_file(params_file_trenches)
    # Skip px_mu conversion
    trench_params = edit_params(trench_params)

    params_segmentation = read_params_file(params_file_segmentation)
    params_segmentation = params_segmentation["mKate2"]["parameters"]

    ###

    # NOTE: Define the GUIs inside the main function so that we can load the params from a file
    # This feels really clunky, defining a function within a function.
    # For one, it prohibits re-use outside of this function.
    # Is there no other way?

    ### Make the Trench GUI

    @magicgui(
        auto_call=False,
        layout="vertical",
        call_button="Calculate",
        tr_height={"min": 0, "max": 100},
        tr_width={"min": 0, "max": 100},
        tr_spacing={
            "min": 0,
            "max": 1000,
        },  # "decimals": 3}, # FIXME how to set the number of decimals in new version of MagicGUI?
        pad_left={"min": 0, "max": 10},
        pad_right={"min": 0, "max": 10},
        pad_top={"min": 0, "max": 20},
        pad_bottom={"min": 0, "max": 20},
    )
    def run_trench_detection_napari(
        channel: Image,
        tr_height: float = trench_params["tr_height"],
        tr_width: int = trench_params["tr_width"],
        tr_spacing: float = trench_params["tr_spacing"],
        pad_left: int = trench_params["pad_left"],
        pad_right: int = trench_params["pad_right"],
        pad_top: int = trench_params["pad_top"],
        pad_bottom: int = trench_params["pad_bottom"],
    ) -> Labels:
        kwargs = {
            # "tr_height"  : tr_height,
            "tr_height": int(tr_height / px_mu),
            "tr_width": tr_width,
            # "tr_width"   : int(tr_width / px_mu),
            "tr_spacing": tr_spacing,
            "pad_left": pad_left,
            "pad_top": pad_top,
            "pad_right": pad_right,
            "pad_bottom": pad_bottom,
        }

        img = channel.data.T

        # find_trenches() returns a multi-dimensional numpy list (and not a numpy array,
        # because each row might have a different number of trenches,
        # after taking into account out-of-bounds issues).
        boxes = generate_boxes(
            kwargs["tr_width"], kwargs["tr_height"], kwargs["tr_spacing"], img.shape[0]
        )

        global trenches
        (trenches, y_offsets) = find_trenches(img, boxes, **kwargs)

        # FIXME or height, width?
        canvas = numpy.zeros(shape=(img.shape[0], img.shape[1]), dtype=numpy.uint16)

        # For now, share the same canvas
        # It would be nice if each row of trenches had its own canvas!
        for r in trenches:
            canvas = regions_arr_to_canvas(canvas, r)

        canvas = canvas.T

        return Labels(canvas, name="Trenches")

    @magicgui(
        auto_call=False,
        layout="vertical",
        call_button="Calculate",
        niblack_k={"min": -1.0, "max": 1.0},
        # Horizontal window size
        # Must be an odd number! TODO How to force or enforce this constraint?
        # Maybe set the value to an odd number, and the step to 2.
        niblack_w_h={"min": 1, "max": 21, "value": 11, "step": 2},
        # Vertical window size
        # Must be an odd number! TODO How to force or enforce this constraint?
        niblack_w_v={"min": 1, "max": 41, "value": 11, "step": 2},
        otsu_multiplier={"max": 10.0},
        garbage_otsu_value={"min": 0, "max": 65536},
        min_size={"min": 0, "max": 1000},
        unsharp_mask_radius={"min": 1, "max": 10},
        unsharp_mask_amount={"min": 1, "max": 10},
    )
    def run_fluor_sharpen_segmentation_napari(
        Layer_Fluor: Image,
        # Layer_Trenches: Labels,
        niblack_k: float = params_segmentation["niblack_k"],
        # NOTE: assumes that the params are passed in as a list of 2 values
        niblack_w_h: int = params_segmentation["niblack_w"][0],
        niblack_w_v: int = params_segmentation["niblack_w"][1],
        otsu_multiplier: float = params_segmentation["otsu_multiplier"],
        garbage_otsu_value: int = params_segmentation["garbage_otsu_value"],
        min_size: float = params_segmentation["min_size"],
        unsharp_mask_radius: int = params_segmentation["unsharp_mask_radius"],
        unsharp_mask_amount: int = params_segmentation["unsharp_mask_amount"],
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

        # FIXME or height, width?
        dim_1 = int(Layer_Fluor.extent.data[1][0]) + 1
        dim_2 = int(Layer_Fluor.extent.data[1][1]) + 1
        canvas = numpy.zeros((dim_2, dim_1), dtype=numpy.uint16)

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

            # Temporarily assign values to any variables which aren't connected to magicgui
            kwargs = {
                "niblack_k": niblack_k,
                "niblack_w": [niblack_w_h, niblack_w_v],
                "otsu_multiplier": otsu_multiplier,
                "garbage_otsu_value": garbage_otsu_value,
                "min_size": min_size,
                "unsharp_mask_radius": unsharp_mask_radius,
                "unsharp_mask_amount": unsharp_mask_amount,
            }

            stack_fl = stack[..., 0]

            # Call the function & get a labeled image
            labeled = fluor_sharpen_segmentation(stack_fl, **kwargs)

            canvas = regions_arr_to_canvas_F_order(canvas, regions, labeled)

        # Done writing to canvas
        canvas = canvas.T
        return Labels(canvas, name="Segmentation mask")

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

        # add_image_layers(images_file, file_nodes, width, height, fields_of_view, frames, z_levels, viewer, channels, layer_params, corrections_file)

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

        # TODO Assume that at this point we will have access to the rows which were returned by the previous widget.

        # Trenches
        # add the gui to the viewer as a dock widget
        viewer.window.add_dock_widget(
            run_trench_detection_napari, name="Trenches", area="left"
        )

        # if a layer gets added or removed, refresh the dropdown choices
        viewer.layers.events.inserted.connect(run_trench_detection_napari.reset_choices)
        viewer.layers.events.removed.connect(run_trench_detection_napari.reset_choices)

        # Segmentation
        # add the gui to the viewer as a dock widget
        viewer.window.add_dock_widget(
            run_fluor_sharpen_segmentation_napari, name="Segmentation", area="left"
        )

        # if a layer gets added or removed, refresh the dropdown choices
        viewer.layers.events.inserted.connect(
            run_fluor_sharpen_segmentation_napari.reset_choices
        )
        viewer.layers.events.removed.connect(
            run_fluor_sharpen_segmentation_napari.reset_choices
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
    default="HDF5/data.h5",
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
