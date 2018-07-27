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
from util import summarize_filenames, get_one
import common
from workflow import get_nd2_frame_anyargs, get_trench_image, get_trench_set_image
import numbers
from cytoolz import get_in, compose

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


def _trench_img(img):
    return RevImage(img.T).options(
        height=120, width=400
    )  # .options(height=10, width=100)#.options(height=int(img.shape[1]*0.6), width=int(img.shape[0]*0.4))


def RevRGB(img, **kwargs):
    return _RevImage(hv.RGB, img, **kwargs)


def show_plot_stack(diags, keys=None):
    if keys is None:
        keys = get_one(diags).keys()
    elif isinstance(keys, str):
        keys = [keys]
    plots = []
    for key in keys:
        plot = hv.HoloMap({t: diag[key] for t, diag in enumerate(diags)}).options(
            title_format=key, finalize_hooks=[_set_active_tool]
        )
        plots.append(plot)
    return hv.Layout.from_values(plots).cols(1).options(normalize=False)


def show_plot_browser(plots, key=None, stream=None, range_xy=None, **kwargs):
    if key is not None:
        if isinstance(key, str):
            key = key.split(".")
        if stream is None:
            plots = get_in(key, plots)
        else:
            plots = compose(partial(get_in, key), plots)
    to_display = {}
    if stream is None:
        initial_plots = plots
    else:
        initial_plots = plots(**stream.contents)
    if range_xy is True:
        range_xy = hv.streams.RangeXY()
    browser = plot_browser(initial_plots, to_display)
    display(browser)
    display_plot_browser_contents(
        plots, to_display, stream=stream, range_xy=range_xy, **kwargs
    )
    return browser


# FROM: https://stackoverflow.com/questions/50415434/how-to-set-active-tools-in-holoviews
def _set_active_tool(plot, element):
    plot.state.toolbar.active_scroll = plot.state.tools[2]


def display_plot_browser_item(
    output, obj, path, stream=None, range_xy=None, regrid=True
):
    if range_xy is None:
        range_xy = hv.streams.RangeXY()
    if isinstance(obj, hv.core.dimension.ViewableElement):
        obj = obj.options(finalize_hooks=[_set_active_tool])
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
                partial(callback, path), streams=[stream, range_xy]
            )  # .collate()
            with output:
                display(dmap)
        else:
            # normal python type
            def callback(o, p, **kwargs):
                with o:
                    clear_output()
                    display(recursive_getattr(plots(**kwargs), p))

            callback(output, path, **stream.contents)
            # we need to do the partial trick, or else output is only bound to the last output of the for loop
            stream.add_subscriber(partial(callback, output, path))


def display_plot_browser_contents(
    plots, to_display, stream=None, range_xy=None, regrid=True
):
    for output, (obj, path) in to_display.items():
        display_plot_browser_item(
            output, obj, path, stream=stream, range_xy=range_xy, regrid=regrid
        )


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
            index_df = df.index.to_frame(index=False)
            row = index_df.iloc[event["new"][0]]  # use first selected row
            params = {}
            for column in stream._df.columns:
                params[column] = row[column]
            stream.event(**params)

        qg.on("selection_changed", handle_selection_changed)

        def update_selection(**kwargs):
            for column in df.index.names:
                if column not in kwargs:
                    return  # stream is missing columns in index
            # idx = df.index.get_loc(tuple(getattr(stream, column) for column in df.index.names))
            # qg._change_selection([idx], 'api', send_msg_to_js=True)
            # the following is equivalent to the above
            idx = tuple(getattr(stream, column) for column in df.index.names)
            qg.change_selection([idx])

        stream.add_subscriber(update_selection)
    # TODO: we can't update selection until widget has already been displayed
    # update_selection(**stream.contents)
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


class DataFrameStream(Stream):
    @classmethod
    def define(cls, name, df):
        params = {"name": param.Parameter(default=name)}
        for column, dtype in df.dtypes.iteritems():
            params[column] = param.Parameter(default=df.iloc[0][column], constant=True)
        params["_df"] = param.Parameter(default=df)
        params["_options"] = {}
        return type(name, (DataFrameStream,), params)

    def __init__(self, **kwargs):
        super(DataFrameStream, self).__init__(**kwargs)
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


class MultiIndexStream(Stream):
    @classmethod
    def define(cls, name, index):
        params = {"name": param.Parameter(default=name)}
        for col_idx, column in enumerate(index.names):
            default = index.levels[col_idx][index.labels[col_idx][0]]
            params[column] = param.Parameter(default=default, constant=True)
        # params['_index'] = param.Parameter(default=index)
        params["_index"] = index
        params["_options"] = {}
        return type(name, (MultiIndexStream,), params)

    def __init__(self, **kwargs):
        super(MultiIndexStream, self).__init__(**kwargs)
        self.transform()  # set self._options

    def transform(self):
        index = self._index
        # TODO
        # mask = np.ones(len(df)).astype(np.bool_)
        index_subset = index
        new_params = {}
        for col_idx, column in enumerate(index.names):
            # TODO
            # options = df[mask][column].unique().tolist()
            options = index_subset._get_level_values(col_idx, unique=True).tolist()
            self._options[column] = options
            value = getattr(self, column)
            if value not in options:
                # TODO
                # value = df[mask][column].iloc[0]
                value = options[0]
                new_params[column] = value
            # TODO
            # mask = (df[column] == value) & mask
            key = (slice(None),) * col_idx + (value,)
            idxs = index_subset.get_locs(key)
            index_subset = index_subset[idxs]
        return new_params

    def __repr__(self):
        cls_name = self.__class__.__name__
        kwargs = ",".join(
            "%s=%r" % (k, v)
            for (k, v) in self.get_param_values()
            if k not in ("name", "_index")
        )
        kwargs += ",_index=<MultiIndex ({} rows)>".format(len(self._index))
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


def dataframe_browser(stream):
    browsers = []
    # TODO
    # for column, dtype in stream._df.dtypes.iteritems():
    #    if issubclass(dtype.type, numbers.Number):
    for column, options in stream._options.items():
        if issubclass(options[0].__class__, numbers.Number):
            kwargs = {"widget": widgets.SelectionSlider}
        else:
            kwargs = {}
        browsers.append(column_browser(column, stream, **kwargs))
    return widgets.VBox(browsers)


def viewer(callback, *streams):
    def callback_wrapper(**kwargs):
        # TODO???
        # kwargs.update({column: kwargs[column] for column in kwargs['_df'].columns if column in kwargs})
        # del kwargs['_df']
        return callback(**kwargs)

    dmap = hv.DynamicMap(callback_wrapper, streams=list(streams))
    return dmap


def dict_viewer(d, *streams, wrapper=None):
    key0 = next(iter(d.keys()))
    fields = key0._fields
    cls = key0.__class__

    def callback(**kwargs):
        key = cls(*[kwargs[field] for field in fields])
        val = d.get(key, None)
        if wrapper:
            val = wrapper(key, val)
        return val

    return viewer(callback, *streams)


def image_viewer(*streams, image_callback=get_nd2_frame_anyargs, regrid=True):
    def callback(x_range, y_range, **kwargs):
        img = image_callback(**kwargs)
        if not isinstance(img, hv.ViewableElement):
            img = RevImage(img)
        if regrid:
            img = datashader.regrid.instance(
                dynamic=False, x_range=x_range, y_range=y_range, aggregator="first"
            )(img).redim.range(z=(0, img.data.max()))
        return img

    return viewer(callback, hv.streams.RangeXY(), *streams)


def trench_viewer(
    trench_bboxes, *streams, channel=None, image_callback=get_trench_image, regrid=False
):
    # TODO: accept trench_bboxes from stream
    def callback(filename, position, t, trench_set, trench, **kwargs):
        # accept channel either from stream or from trench_viewer kwargs
        _channel = kwargs.get("channel", None) or channel
        return _trench_img(
            image_callback(
                trench_bboxes, filename, position, _channel, t, trench_set, trench
            )
        )

    return image_viewer(*streams, image_callback=callback, regrid=regrid)


def trench_set_viewer(
    trench_bboxes,
    *streams,
    channel=None,
    image_callback=get_trench_set_image,
    regrid=True,
):
    # TODO: accept trench_bboxes from stream
    def callback(filename, position, t, trench_set, trench, **kwargs):
        # accept channel either from stream or from trench_viewer kwargs
        _channel = kwargs.get("channel", None) or channel
        return RevImage(
            image_callback(trench_bboxes, filename, position, _channel, t, trench_set)
        )

    return image_viewer(*streams, image_callback=callback, regrid=regrid)


def show_frame_info(df, stream):
    df_noindex = df.reset_index()
    qg = show_grid(
        df_noindex.iloc[0],
        header_height=20,
        defaultColumnWidth=300,
        forceFitColumns=True,
    )

    def update_frame_info(**kwargs):
        idx = df.index.get_loc(
            tuple(getattr(stream, column) for column in df.index.names)
        )
        new_df = df_noindex.iloc[idx]
        qg._df.loc[:, 0] = new_df
        qg._unfiltered_df.loc[:, 0] = new_df
        qg._update_table(triggered_by="cell_change", fire_data_change_event=True)

    stream.add_subscriber(update_frame_info)
    update_frame_info(**stream.contents)
    return qg
