import numpy as np
import re
from bokeh.layouts import column, gridplot
from bokeh.plotting import figure, show
from bokeh.models import Range1d, LinearAxis, BoxAnnotation
from bokeh.models.callbacks import CustomJS
from functional import compose
from align import align_mafft, read_ab1, trim_to_ref
from util import grouper


def show_chromatogram_alignment(
    ab1_files,
    ref_seq,
    aligner=compose(trim_to_ref, align_mafft),
    margin=10,
    detail_width=100,
    plot_width=900,
    base_width=10,
    detail_base_width=1,
    text_size=12,
    text_margin=2,
    row_margin=5,
):
    text_height = text_size + text_margin
    ab1s = [read_ab1(ab1_file) for ab1_file in ab1_files]
    sequences = {"ref": ref_seq, **{ab1["name"]: ab1["sequence"] for ab1 in ab1s}}
    msa = aligner(sequences)
    # calculate dimensions
    sequence_strs = [str(m.seq) for m in msa]
    max_length = max(len(seq) for seq in sequence_strs)
    bases_per_row = plot_width // base_width
    grouped_sequences = [
        [np.array(s) for s in grouper(seq, bases_per_row)] for seq in sequence_strs
    ]
    num_rows = max(len(row) for row in grouped_sequences)
    num_seqs = len(sequence_strs)
    row_height = text_height * num_seqs + row_margin
    # make figures
    overview_boxes = []
    overview_x_range = Range1d(
        -base_width, bases_per_row + base_width
    )  # is margin right?
    overview_y_range = Range1d(row_height * num_rows, 0)
    overview_plot = figure(
        plot_width=plot_width,
        plot_height=row_height * num_rows,
        tools="tap,save",
        x_range=overview_x_range,
        y_range=overview_y_range,
        output_backend="webgl",
    )
    overview_plot.axis.visible = False
    overview_plot.min_border = 0
    overview_plot.min_border_bottom = 0
    overview_plot.grid.visible = False
    overview_plot.outline_line_color = None
    for row in range(num_rows):
        for seq in range(num_seqs):
            sequence = grouped_sequences[seq][row]
            r = overview_plot.text(
                x=np.arange(len(sequence)) + 0.5,
                y=text_height * (seq + 1) + row_height * row,
                text=[s for s in sequence],
                text_align="center",
                text_baseline="bottom",
                text_font="Courier",
                text_font_size="{}px".format(text_size),
                text_color=np.where(
                    np.logical_and(
                        sequence != grouped_sequences[0][row], sequence != "-"
                    ),
                    "red",
                    "black",
                ),
            )
            if seq == 0:
                r.glyph.text_font_style = "bold"
        box = BoxAnnotation(
            left=0,
            right=0,
            top=row_height * row,
            bottom=row_height * row + text_height * num_seqs,
            fill_color="#ffd3f1",
            fill_alpha=0.5,
        )
        box.visible = False
        overview_plot.add_layout(box)
        overview_boxes.append(box)
    detail_plots = []
    x_min = min(np.min(ab1["peak_calls"]) for ab1 in ab1s)
    x_max = max(np.max(ab1["peak_calls"]) for ab1 in ab1s)
    detail_x_range = Range1d(
        x_min,
        x_min + detail_width,
        bounds=(x_min - margin, x_max + margin),
        max_interval=detail_width,
    )
    y_margin_factor = 0.6
    detail_y_range = Range1d(-y_margin_factor, 1)
    plot_ref_sequence = figure(
        plot_width=plot_width,
        plot_height=int(50 // (1 + y_margin_factor)),
        output_backend="webgl",
        # tools=['xpan', 'reset','xwheel_zoom', BoxZoomTool(dimensions='width')],
        tools=["xpan"],
        active_drag="xpan",
        # active_scroll='xwheel_zoom',
        # title=name,
        x_range=detail_x_range,
        y_range=Range1d(-y_margin_factor, 0),
    )
    r = plot_ref_sequence.text(
        x=range(len(sequence_strs[0])),
        y=-0.3,
        text=[s for s in sequence_strs[0]],
        text_align="center",
        text_baseline="middle",
        text_font_size="16px",
    )
    r.glyph.text_font_style = "bold"
    plot_ref_sequence.min_border = 0
    plot_ref_sequence.grid.visible = False
    plot_ref_sequence.outline_line_color = None
    detail_plots.append(plot_ref_sequence)
    for aligned_sequence, ab1 in zip(sequence_strs[1:], ab1s):
        name, sequence, peaks, peak_calls, quality = map(
            ab1.__getitem__, ("name", "sequence", "peaks", "peak_calls", "quality")
        )
        peak_normalization = max(p.max() for p in peaks.values())
        p = figure(
            plot_width=plot_width,
            plot_height=50,
            output_backend="webgl",
            # tools=['xpan', 'reset','xwheel_zoom', BoxZoomTool(dimensions='width')],
            tools=["xpan"],
            active_drag="xpan",
            # active_scroll='xwheel_zoom',
            # title=name,
            x_range=detail_x_range,
            y_range=detail_y_range,
        )
        p.extra_y_ranges = {"q": Range1d(start=-y_margin_factor * 100, end=100)}
        p.add_layout(LinearAxis(y_range_name="q"), "right")
        local_start = 0
        for match in re.finditer(r"([^-]+)", aligned_sequence):
            aligned_start = match.start()
            aligned_end = match.end()
            sequence_section = aligned_sequence[aligned_start:aligned_end]
            section_length = len(sequence_section)
            local_end = local_start + section_length - 1
            quality_section = quality[local_start : local_end + 1]
            quality_x = np.linspace(
                aligned_start, aligned_end, len(quality_section), endpoint=False
            )
            p.rect(
                x=quality_x,
                y=quality_section / 2,
                width=1,
                height=quality_section,
                y_range_name="q",
                color="#e7ecf6",
            )
            for color, base in zip(
                ["red", "green", "blue", "black"], ["T", "A", "C", "G"]
            ):
                peak = peaks[base]
                peak_left_start = int(
                    (peak_calls[local_start - 1] + peak_calls[local_start]) / 2
                )
                peak_left = peak[peak_left_start : peak_calls[local_start]]
                peak_left_x = np.linspace(
                    aligned_start - 0.5, aligned_start, len(peak_left), endpoint=False
                )
                peak_section = peak[peak_calls[local_start] : peak_calls[local_end]]
                peak_section_x = np.linspace(
                    aligned_start, aligned_end - 1, len(peak_section), endpoint=False
                )
                if local_end + 1 < len(peak_calls):
                    peak_right_end = int(
                        (peak_calls[local_end] + peak_calls[local_end + 1]) / 2
                    )
                else:
                    peak_right_end = peak_calls[local_end]
                peak_right = peak[peak_calls[local_end] : peak_right_end]
                peak_right_x = np.linspace(
                    aligned_end - 1, aligned_end - 0.5, len(peak_right), endpoint=False
                )
                peak_x = np.concatenate((peak_left_x, peak_section_x, peak_right_x))
                peak_y = np.concatenate((peak_left, peak_section, peak_right))
                p.line(peak_x, peak_y / peak_normalization, line_color=color)
            local_start += len(sequence_section)
        detail_plots.append(p)
        # p2 = figure(plot_width=plot_width, plot_height=30, tools=['xpan'], active_drag='xpan', x_range=detail_x_range,
        #             output_backend='webgl')
        aligned_sequence_ary = np.array([s for s in aligned_sequence])
        r = p.text(
            x=range(len(aligned_sequence)),
            y=-0.3,
            text=[s for s in aligned_sequence],
            text_align="center",
            text_baseline="middle",
            text_font_size="16px",
            text_color=np.where(
                np.logical_and(
                    aligned_sequence_ary != np.array([s for s in sequence_strs[0]]),
                    aligned_sequence_ary != "-",
                ),
                "red",
                "black",
            ),
        )
        # detail_plots.append(p2)
        p.min_border = 0
        # p2.min_border = 0
        p.grid.visible = False
        # p2.grid.visible = False
        p.outline_line_color = None
        # p2.outline_line_color = None
    for p in detail_plots:
        p.xaxis.visible = False
    for p in detail_plots:
        p.yaxis.visible = False
    detail_plots[-1].height = 80
    detail_plots[-1].xaxis.visible = True
    detail_plots[-1].outline_line_color = None
    detail_plots[-1].xaxis.axis_line_color = None
    # callback
    BOX_JS = """
    for (i=0; i<boxes.length; i++) {
        var start_in_row = false;
        var end_in_row = false;
        if ((i*bases_per_row <= start) && (start < (i+1)*bases_per_row)) {
           start_in_row = true;
        }
        if ((i*bases_per_row <= end) && (end < (i+1)*bases_per_row)) {
           end_in_row = true;
        }
        //console.log('row',i,'startinrow',start_in_row,'endinrow',end_in_row);
        if (start_in_row || end_in_row) {
           boxes[i].visible = true;
        }
        else {
           boxes[i].visible = false;
        }
        if (start_in_row) {
           boxes[i].left = (start-i*bases_per_row);
           boxes[i].right = bases_per_row;
        }
        else {
           boxes[i].left = 0;
        }
        if (end_in_row) {
           boxes[i].left = 0;
           boxes[i].right = (end-i*bases_per_row);
        }
        else {
           boxes[i].right = bases_per_row;
        }
        if ((i*bases_per_row > start) && (end >= (i+1)*bases_per_row)) {
           boxes[i].visible = true;
           boxes[i].left = 0;
           boxes[i].right = bases_per_row;
        }
    }
    """
    detail_scroll_callback = CustomJS(
        args=dict(boxes=overview_boxes, x_range=detail_plots[0].x_range),
        code="""
                                           var bases_per_row = {bases_per_row};
                                           var start = x_range.start;
                                           var end = x_range.end;
                                           """.format(
            bases_per_row=bases_per_row
        )
        + BOX_JS,
    )
    # TODO!!!
    detail_plots[0].x_range.js_on_change("end", detail_scroll_callback)
    range_callback = CustomJS(
        args=dict(boxes=overview_boxes, plot=detail_plots[0]),
        code="""
                                   var bases_per_row = {bases_per_row};
                                   var detail_width = {detail_width};
                                   var detail_base_width = {detail_base_width};
                                   var length = {length};
                                   console.log('x', cb_obj.x, 'y', cb_obj.y);
                                   var row = Math.floor(cb_obj.y / {row_height});
                                   var start = Math.max(Math.floor(cb_obj.x) + row*bases_per_row - Math.floor(detail_width / 2), 0);
                                   var end = start + detail_width; // also bound
                                   if (end > length) {{
                                       end = length;
                                       start = end - detail_width;
                                   }}
                                   //plot.trigger('change');
                                   plot.x_range.start = start;
                                   plot.x_range.end = end;
                                   """.format(
            bases_per_row=bases_per_row,
            detail_width=detail_width,
            length=max_length,
            row_height=row_height,
            detail_base_width=detail_base_width,
        ),
    )
    overview_plot.js_on_event("tap", range_callback)
    # initialize boxes
    end = int(detail_plots[0].x_range.end)
    for i, box in enumerate(overview_boxes):
        box.visible = True
        if (i * bases_per_row <= end) and (end < (i + 1) * bases_per_row):
            box.left = 0
            box.right = end - i * bases_per_row
            break
        else:
            box.left = 0
            box.right = bases_per_row
    return column(overview_plot, gridplot([[p] for p in detail_plots]))


def show_unrolled_chromatogram(
    ab1_files,
    ref_seq,
    aligner=compose(trim_to_ref, align_mafft),
    margin=10,
    detail_width=100,
    plot_width=900,
    base_width=10,
    detail_base_width=10,
    text_size=12,
    text_margin=2,
    row_margin=5,
):
    text_height = text_size + text_margin
    ab1s = [read_ab1(ab1_file) for ab1_file in ab1_files]
    sequences = {"ref": ref_seq, **{ab1["name"]: ab1["sequence"] for ab1 in ab1s}}
    msa = aligner(sequences)
    # calculate dimensions
    sequence_strs = [np.array(str(m.seq)) for m in msa]
    max_length = max(len(seq) for seq in sequence_strs)
    bases_per_row = plot_width // base_width
    grouped_sequences = [list(grouper(seq, bases_per_row)) for seq in sequence_strs]
    num_rows = max(len(row) for row in grouped_sequences)
    num_seqs = len(sequence_strs)
    row_height = text_height * num_seqs + row_margin
    # make figures
    # overview_plots[-1].height += 30
    for row in range(num_rows):
        for seq in range(num_seqs):
            pass
    detail_plots = []
    detail_x_range = None
    for aligned_sequence, ab1 in zip(sequence_strs[1:], ab1s):
        name, sequence, peaks, peak_calls, quality = map(
            ab1.__getitem__, ("name", "sequence", "peaks", "peak_calls", "quality")
        )
        if detail_x_range is None:
            x_min = np.min(peak_calls)
            x_max = np.max(peak_calls)
            detail_x_range = Range1d(
                x_min,
                x_min + detail_width,
                bounds=(x_min - margin, x_max + margin),
                max_interval=detail_width,
            )
        p = figure(
            plot_width=plot_width,
            plot_height=50,
            output_backend="webgl",
            # tools=['xpan', 'reset','xwheel_zoom', BoxZoomTool(dimensions='width')],
            tools=["xpan"],
            active_drag="xpan",
            # active_scroll='xwheel_zoom',
            # title=name,
            x_range=detail_x_range,
            # y_range=y_range
        )
        p.extra_y_ranges = {"q": Range1d(start=0, end=100)}
        p.add_layout(LinearAxis(y_range_name="q"), "right")
        local_start = 0
        for match in re.finditer(r"([^-]+)", aligned_sequence):
            aligned_start = match.start()
            aligned_end = match.end()
            sequence_section = aligned_sequence[aligned_start:aligned_end]
            section_length = len(sequence_section)
            local_end = local_start + section_length - 1
            quality_section = quality[local_start : local_end + 1]
            # print('qs',len(quality_section),'adist',aligned_start,aligned_end-1,'diff',aligned_end-1-aligned_start)
            quality_x = np.linspace(
                aligned_start, aligned_end, len(quality_section), endpoint=False
            )
            p.rect(
                x=quality_x,
                y=quality_section / 2,
                width=1,
                height=quality_section,
                y_range_name="q",
                color="#e7ecf6",
            )
            for color, base in zip(
                ["red", "green", "blue", "black"], ["T", "A", "C", "G"]
            ):
                peak = peaks[base]
                # peak_section = peak[peak_calls[local_start]:peak_calls[local_end]]
                # p.line(np.linspace(aligned_start, aligned_end-1, len(peak_section)), peak_section, line_color=color)
                peak_left_start = int(
                    (peak_calls[local_start - 1] + peak_calls[local_start]) / 2
                )
                peak_left = peak[peak_left_start : peak_calls[local_start]]
                peak_left_x = np.linspace(
                    aligned_start - 0.5, aligned_start, len(peak_left), endpoint=False
                )
                peak_section = peak[peak_calls[local_start] : peak_calls[local_end]]
                peak_section_x = np.linspace(
                    aligned_start, aligned_end - 1, len(peak_section), endpoint=False
                )
                if local_end + 1 < len(peak_calls):
                    peak_right_end = int(
                        (peak_calls[local_end] + peak_calls[local_end + 1]) / 2
                    )
                else:
                    peak_right_end = peak_calls[local_end]
                peak_right = peak[peak_calls[local_end] : peak_right_end]
                peak_right_x = np.linspace(
                    aligned_end - 1, aligned_end - 0.5, len(peak_right), endpoint=False
                )
                # p.line(peak_section_x, peak_section, line_color=color)
                peak_x = np.concatenate((peak_left_x, peak_section_x, peak_right_x))
                peak_y = np.concatenate((peak_left, peak_section, peak_right))
                p.line(peak_x, peak_y, line_color=color)
            local_start += len(sequence_section)
        detail_plots.append(p)
        p2 = figure(
            plot_width=plot_width,
            plot_height=30,
            tools=["xpan"],
            active_drag="xpan",
            x_range=detail_x_range,
            output_backend="webgl",
        )
        r = p2.text(
            x=range(len(aligned_sequence)),
            y=0.1,
            text=[s for s in aligned_sequence],
            text_align="center",
            text_baseline="middle",
        )
        # r.glyph.text_font_style='bold'
        detail_plots.append(p2)
        p.min_border = 0
        p2.min_border = 0
        p.grid.visible = False
        p2.grid.visible = False
        p.outline_line_color = None
        p2.outline_line_color = None
    for p in detail_plots:
        p.xaxis.visible = False
    for p in detail_plots:
        p.yaxis.visible = False
    detail_plots[-1].height = 50
    detail_plots[-1].xaxis.visible = True
    show(plot)
