import numpy as np
import holoviews as hv
import ipywidgets as widgets

# from holoviews.operation.datashader import aggregate, datashade, dynspread, shade, regrid
from holoviews.operation import datashader
from IPython.display import display, clear_output, HTML
from holoviews.streams import Stream, param
from matplotlib.colors import hex2color
import qgrid
from collections.abc import Mapping, Sequence
from functools import partial, reduce
import uuid
from util import iterate_getattr, summarize_filenames
import common
from workflow import get_nd2_frame

# TODO
channel_to_color = {
    "BF": "#ffffff",
    "MCHERRY": "#e22400",
    "GFP": "#76ba40",
    "CY5": "#e292fe",
    "BFP": "#3a87fd",
    "YFP": "#f5eb00",
}


def RevImage(img, **kwargs):
    return _RevImage(hv.Image, img, **kwargs)


def _RevImage(cls, img, **kwargs):
    return cls(img[::-1], bounds=(0, 0, img.shape[1], img.shape[0])).options(
        invert_yaxis=True
    )


def RevRGB(img, **kwargs):
    return _RevImage(hv.RGB, img, **kwargs)


def display_plot_browser(plots, stream=None, **kwargs):
    to_display = {}
    if stream is None:
        initial_plots = plots
    else:
        initial_plots = plots(**stream.contents)
    browser = plot_browser(initial_plots, to_display)
    display(browser)
    display_plot_browser_contents(plots, to_display, stream=stream, **kwargs)
    return browser


def display_plot_browser_contents(plots, to_display, stream=None, regrid=True):
    for output, (obj, path) in to_display.items():
        if stream is None:
            with output:
                if isinstance(obj, hv.core.dimension.ViewableElement):
                    if regrid:
                        obj = obj.map(
                            lambda img: datashader.regrid(
                                img, aggregator="first"
                            ).redim.range(z=(0, img.data.max())),
                            lambda obj: isinstance(obj, (hv.Image, hv.RGB, hv.Raster)),
                        )
                        if hasattr(obj, "collate"):
                            obj = obj.collate()
                display(obj)
        else:
            # print(path, obj.__class__)
            if isinstance(
                obj, hv.core.dimension.ViewableElement
            ):  # TODO: is this the right comparison?

                def callback(p, x_range, y_range, **kwargs):
                    plot = recursive_getattr(plots(**kwargs), p)
                    # allow customizable aggregators
                    return plot.map(
                        lambda img: datashader.regrid.instance(
                            dynamic=False,
                            x_range=x_range,
                            y_range=y_range,
                            aggregator="first",
                        )(img).redim.range(z=(0, img.data.max())),
                        lambda obj: isinstance(obj, (hv.Image, hv.RGB, hv.Raster)),
                    )

                dmap = hv.DynamicMap(
                    partial(callback, path), streams=[stream, hv.streams.RangeXY()]
                )  # .collate()
                with output:
                    display(dmap)
            else:
                # normal python type
                def callback(o, p, **kwargs):
                    with o:
                        clear_output()
                        display(recursive_getattr(plots(**kwargs), p))
                        # display((p, kwargs, np.random.random(), recursive_getattr(plots(**kwargs), p)))
                        # display((p, np.random.random()))

                callback(output, path, **stream.contents)
                # we need to do the partial trick, or else output is only bound to the last output of the for loop
                stream.add_subscriber(partial(callback, output, path))


def plot_browser(plots, to_display=None, path=()):
    if not isinstance(plots, Mapping):
        raise NotImplementedError
    children = []
    singleton_children = [k for k, v in plots.items() if not isinstance(v, Mapping)]
    for k in singleton_children:
        label = widgets.HTML("<b>{}</b>".format(k))
        output = widgets.Output()
        to_display[output] = (plots[k], path + (k,))
        child = widgets.HBox([label, output])
        children.append(child)
    nested_children = [k for k in plots.keys() if k not in singleton_children]
    if nested_children:
        accordion_children = [
            plot_browser(plots[k], to_display=to_display, path=(path + (k,)))
            for k in nested_children
        ]
        accordion = widgets.Accordion(children=accordion_children)
        for i, k in enumerate(nested_children):
            accordion.set_title(i, k)
        children.append(accordion)
    return widgets.VBox(children)


def resize_qgrid_header(widget, height):
    widget_class = "qgrid-{}".format(uuid.uuid1())
    widget.add_class(widget_class)
    filter_height = height + 12  # TODO: why 12??
    css = """
    <style>
    .{widget_class} .slick-header-column.ui-state-default {{
        height: {height}px;
    }}
    .{widget_class} .slick-column-name {{
        width: {height}px;
        overflow-wrap: break-word;
    }}
    .{widget_class} .slick-column-name {{
        white-space: normal;
    }}
    .{widget_class} .filter-button {{
        height: {filter_height}px;
    }}
    </style>"""
    display(
        HTML(
            css.format(
                widget_class=widget_class, height=height, filter_height=filter_height
            )
        )
    )


Selected = Stream.define("Selected", selected=None)


def show_grid(df, header_height=100, stream=None, **kwargs):
    qg = qgrid.show_grid(
        df,
        grid_options={
            "forceFitColumns": False,
            "editable": False,
            "enableColumnReorder": True,
            "defaultColumnWidth": 90,
            **kwargs,
        },
        precision=1,
    )
    resize_qgrid_header(qg, header_height)
    if stream is not None:

        def handle_selection_changed(event, widget):
            stream.event(selected=df.index[event["new"][0]])

        qg.on("selection_changed", handle_selection_changed)
        # TODO: update qgrid selection on stream event (without recursion!)
    return qg


def _select(xs, mask):
    if len(xs) != len(mask):
        raise ValueError("mask length does not match")
    return [xs[i] for i in range(len(xs)) if mask[i]]


def composite_channels(imgs, hexcolors, scale=True):
    colors = [hex2color(hexcolor) for hexcolor in hexcolors]
    return _composite_channels(imgs, colors, scale=scale)


def _composite_channels(channel_imgs, colors, scale=True):
    if len(channel_imgs) != len(colors):
        raise ValueError("expecting equal numbers of channels and colors")
    num_channels = len(channel_imgs)
    if scale:
        scaled_imgs = [
            channel_imgs[i][:, :, np.newaxis] / np.percentile(channel_imgs[i], 99.9)
            for i in range(num_channels)
        ]
        for scaled_img in scaled_imgs:
            np.clip(scaled_img, 0, 1, scaled_img)  # clip in place
    else:
        scaled_imgs = channel_imgs
    imgs_to_combine = [
        scaled_imgs[i] * np.array(colors[i]) for i in range(num_channels)
    ]
    if not len(imgs_to_combine):
        imgs_to_combine = [np.ones(colored_imgs[0].shape)]  # white placeholder
    img = imgs_to_combine[0]
    for img2 in imgs_to_combine[1:]:
        img = 1 - (1 - img) * (1 - img2)
    return img


def multichannel_selector(frames):
    channels = frames.attrs["metadata"]["channels"]
    num_channels = len(channels)
    # colors = [hex2color(channel_colors[channel]) for channel in channels]
    channel_boxes = []
    channel_widgets = []
    channel_enabled = [True] * num_channels
    channel_colors = [channel_to_color[channel] for channel in channels]
    for i, channel in enumerate(channels):
        solo_button = widgets.Button(
            description="S", layout=widgets.Layout(width="10%")
        )
        enabled_button = widgets.ToggleButton(
            description=channel, value=channel_enabled[i]
        )
        solo_button._button_to_enable = enabled_button
        color_picker = widgets.ColorPicker(concise=True, value=channel_colors[i])
        channel_box = widgets.HBox([solo_button, enabled_button, color_picker])
        channel_widgets.append([solo_button, enabled_button, color_picker, channel_box])
    solo_buttons, enabled_buttons, color_pickers, channel_boxes = zip(*channel_widgets)
    channels_box = widgets.VBox(channel_boxes)
    display_settings_stream = DisplaySettings()

    def update_enabled_channels(change):
        channel_enabled = [button.value for button in enabled_buttons]
        display_settings_stream.event(channel_enabled=channel_enabled)

    def update_solo(solo_button):
        if (
            solo_button._button_to_enable.value
            and sum([b.value for b in enabled_buttons]) == 1
        ):
            for enabled_button in enabled_buttons:
                enabled_button.value = True
        else:
            for enabled_button in enabled_buttons:
                enabled_button.value = enabled_button == solo_button._button_to_enable
        # update_enabled_channels(None)

    for solo_button in solo_buttons:
        solo_button.on_click(update_solo)
    for enabled_button in enabled_buttons:
        enabled_button.observe(update_enabled_channels, names="value")

    def update_channel_colors(change):
        channel_colors = [color_picker.value for color_picker in color_pickers]
        display_settings_stream.event(channel_colors=channel_colors)

    for color_picker in color_pickers:
        color_picker.observe(update_channel_colors, names="value")
    return channels_box, display_settings_stream


# FrameStream = Stream.define('Frame', t=0, v=0)

# def positionpoints_browser(frames, frame_stream):
#     num_positionpoints = len(frames.attrs['metadata']['frames'])
#     #play_buttons = widgets.Play(interval=10, min=0, max=num_positionpoints, step=1)
#     back_step_button = widgets.Button(description='<', layout=widgets.Layout(width='10%'))
#     forward_step_button = widgets.Button(description='>', layout=widgets.Layout(width='10%'))
#     t_slider = widgets.IntSlider(label='t', min=0, max=num_positionpoints, step=1, value=0, continuous_update=False)
#     slider_box = widgets.HBox([back_step_button, t_slider, forward_step_button])
#     t_slider.observe(lambda change: frame_stream.event(t=change['new']), names='value')
#     return slider_box

# def frame_browser(frames, frame_stream):
#     num_positionpoints = len(frames.attrs['metadata']['frames'])
#     num_fovs = len(frames.attrs['metadata']['fields_of_view'])
#     t_slider = widgets.IntSlider(label='t', min=0, max=num_positionpoints, step=1, value=0, continuous_update=False)
#     v_slider = widgets.IntSlider(label='v', min=0, max=num_fovs, step=1, value=0, continuous_update=False)
#     slider_box = widgets.VBox([v_slider, t_slider])
#     t_slider.observe(lambda change: frame_stream.event(t=change['new']), names='value')
#     v_slider.observe(lambda change: frame_stream.event(v=change['new']), names='value')
#     return slider_box


def big_image_viewer(positions, frame_stream=None):
    num_channels = positions[0].shape[0]  # TODO
    if frame_stream is None:
        frame_stream = FrameStream()
    slider_box = frame_browser(positions, frame_stream)
    channels_box, display_settings_stream = multichannel_selector(positions)

    def image_callback(t, v, channel_enabled, channel_colors):
        pos_stack = positions[str(v)]
        channel_imgs = [
            pos_stack[c, t, :, :] for c in range(num_channels) if channel_enabled[c]
        ]
        img = composite_channels(channel_imgs, _select(channel_colors, channel_enabled))
        viewer = RevRGB(img)
        return viewer

    image = hv.DynamicMap(
        image_callback, streams=[frame_stream, display_settings_stream]
    )
    image = datashader.regrid(image)
    image = image.opts(plot={"width": 500, "height": 500})
    output = widgets.Output()
    box = widgets.VBox([widgets.HBox([channels_box, slider_box]), output])
    display(box)
    with output:
        display(image)
    return box


FrameChannels = Stream.define(
    "FrameChannels", channel_enabled=None, channel_colors=None
)


class DataframeStream(Stream):
    @classmethod
    def define(cls, name, df):
        params = {"name": param.Parameter(default=name)}
        for column, dtype in df.dtypes.iteritems():
            params[column] = param.Parameter(default=df.iloc[0][column], constant=True)
        params["_df"] = param.Parameter(default=df)
        params["_options"] = {}
        return type(name, (DataframeStream,), params)

    def __init__(self, **kwargs):
        super(DataframeStream, self).__init__(**kwargs)
        self.transform()  # set self._options

    def transform(self):
        df = self._df
        mask = np.ones(len(df)).astype(np.bool_)
        new_params = {}
        for column in df.columns:
            options = df[mask][column].unique().tolist()
            self._options[column] = options
            value = getattr(self, column)
            if value not in options:
                value = df[mask][column].iloc[0]
                new_params[column] = value
            mask = (df[column] == value) & mask
        return new_params

    def __repr__(self):
        cls_name = self.__class__.__name__
        kwargs = ",".join(
            "%s=%r" % (k, v)
            for (k, v) in self.get_param_values()
            if k not in ("name", "_df")
        )
        kwargs += ",_df=<DataFrame ({} rows)>".format(len(self._df))
        if not self._rename:
            return "%s(%s)" % (cls_name, kwargs)
        else:
            return "%s(%r, %s)" % (cls_name, self._rename, kwargs)


def column_browser(
    name, stream, widget=widgets.Dropdown, format_function=lambda x: map(str, x)
):
    left_arrow = widgets.Button(description="<", layout=widgets.Layout(width="3%"))
    right_arrow = widgets.Button(description=">", layout=widgets.Layout(width="3%"))
    selector = widget(description=name, options=range(1))
    browser = widgets.HBox([left_arrow, selector, right_arrow])

    def increment(inc, button):
        options = stream._options[name]
        idx = options.index(getattr(stream, name))
        new_idx = (idx + inc) % len(options)
        stream.event(**{name: options[new_idx]})

    left_arrow.on_click(partial(increment, -1))
    right_arrow.on_click(partial(increment, 1))

    def set_selection(change):
        if change["name"] == "value":
            # TODO: do I need to do this check?
            if change["new"] != getattr(stream, name):
                stream.event(**{name: change["new"]})

    selector.observe(set_selection)

    def update_selector(**kwargs):
        options = stream._options[name]
        options_for_display = dict(zip(format_function(options), options))
        try:
            selector.options = options_for_display
        except:
            pass
        selector.value = kwargs[name]

    stream.add_subscriber(update_selector)
    return browser
    return position_browser


def frame_browser(df, frame_stream):
    filename_browser = widgets.HBox([left_arrow, dropdown, filename_right])
    box = widgets.VBox([])


def image_viewer(*streams, image_callback=get_nd2_frame):
    def callback(x_range, y_range, **kwargs):
        img = RevImage(image_callback(**kwargs))
        return datashader.regrid.instance(
            dynamic=False, x_range=x_range, y_range=y_range, aggregator="first"
        )(img).redim.range(z=(0, img.data.max()))

    dmap = hv.DynamicMap(callback, streams=streams + [hv.streams.RangeXY()])
    return dmap
