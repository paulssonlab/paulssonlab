import gdstk
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextToPath

from paulssonlab.microfluidics_design.geometry import (
    Cell,
    from_matplotlib_path,
    from_matplotlib_path_codes,
)
from paulssonlab.microfluidics_design.util import memoize


@memoize
def glyph(glyph_id, font_prop, verts, codes, tolerance=0.1, layer=0, datatype=0):
    cell = Cell(f"Glyph {glyph_id}")
    cell.add(
        *from_matplotlib_path_codes(
            verts, codes, tolerance=tolerance, layer=layer, datatype=datatype
        )
    )
    return cell


# FROM: https://github.com/matplotlib/matplotlib/blob/main/lib/matplotlib/textpath.py
def _text_line(s, font_prop, position=(0, 0), ismath=False, layer=0, datatype=0):
    """Convert text *s* to GDS cell references and polygons. Each glyph is
    rendered as a reference so its polygons are only rendered once.

    Parameters
    ----------
    s : str
        The text to be converted.
    font_prop : `matplotlib.font_manager.FontProperties`
        The font properties for the text.
    ismath : {False, True, "TeX"}
        If True, use mathtext parser.  If "TeX", use tex for rendering.

    Returns
    -------
    cells : list
        A list of cell references for each glyph.
    """
    position = np.array(position)
    texttopath = TextToPath()
    font_prop_unscaled = font_prop.copy()
    font_prop_unscaled.set_size(texttopath.FONT_SCALE)
    if ismath == "TeX":
        glyph_info, glyph_map, rects = texttopath.get_glyphs_tex(font_prop_unscaled, s)
    elif not ismath:
        font = texttopath._get_font(font_prop_unscaled)
        glyph_info, glyph_map, rects = texttopath.get_glyphs_with_font(font, s)
    else:
        glyph_info, glyph_map, rects = texttopath.get_glyphs_mathtext(
            font_prop_unscaled, s
        )
    matplotlib_scale = font_prop.get_size_in_points() / texttopath.FONT_SCALE
    glyph_cells = {}
    polys = []
    for glyph_id, x, y, scale in glyph_info:
        offset = np.array([x, y])
        if glyph_id in glyph_cells:
            glyph_cell = glyph_cells[glyph_id]
        else:
            glyph_cell = glyph(
                glyph_id,
                font_prop_unscaled,
                *glyph_map[glyph_id],
                layer=layer,
                datatype=datatype,
            )
            glyph_cells[glyph_id] = glyph_cell
        polys.append(
            gdstk.Reference(
                glyph_cell,
                position + offset * matplotlib_scale,
                magnification=scale * matplotlib_scale,
            )
        )
    for verts1, codes1 in rects:
        polys.add(*from_matplotlib_path_codes(verts1, codes1))
    return polys


# SEE _get_layout at https://github.com/matplotlib/matplotlib/blob/main/lib/matplotlib/text.py
# did not implement special handling of multiline TeX
def _text(s, font_prop, position=(0, 0), ismath=False, layer=0, datatype=0):
    x, y = position
    texttopath = TextToPath()
    # use the characters "lp" to measure descent, so that descent is fixed
    # for all lines
    width, height, descent = texttopath.get_text_width_height_descent(
        "lp", font_prop, False
    )
    line_height = height + descent
    polys = []
    for line_num, line in enumerate(s.split("\n")):
        polys.extend(
            _text_line(
                line,
                font_prop,
                position=(x, y - line_num * height),
                ismath=ismath,
                layer=layer,
                datatype=datatype,
            )
        )
    return polys


def text(
    s,
    size,
    position=(0, 0),
    ismath=False,
    layer=0,
    datatype=0,
):
    font_prop = FontProperties(
        family="DejaVu Sans",
        style="normal",
        variant="normal",
        weight="normal",
        size=size,
        math_fontfamily="dejavusans",
    )
    return _text(
        s, font_prop, position=position, ismath=ismath, layer=layer, datatype=datatype
    )
