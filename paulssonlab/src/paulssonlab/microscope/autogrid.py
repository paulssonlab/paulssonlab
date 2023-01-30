import xml.etree.ElementTree as ET
import click
import numpy as np


def parse_position(e):
    return dict(
        name=e.find("strName").attrib["value"],
        checked=(e.find("bChecked").attrib["value"].lower() != "false"),
        x=float(e.find("dXPosition").attrib["value"]),
        y=float(e.find("dYPosition").attrib["value"]),
        z=float(e.find("dZPosition").attrib["value"]),
        pfs_offset=float(e.find("dPFSOffset").attrib["value"]),
    )


def position_element(idx, name, x, y, z, pfs_offset=-1, checked=True):
    tag = f"Point{idx:05d}"
    position = ET.Element(tag, runtype="NDSetupMultipointListItem")
    ET.SubElement(position, "bChecked", runtype="bool", value=str(checked).lower())
    ET.SubElement(position, "strName", runtype="CLxStringW", value=str(name))
    ET.SubElement(position, "dXPosition", runtype="double", value=f"{x:.15f}")
    ET.SubElement(position, "dYPosition", runtype="double", value=f"{y:.15f}")
    ET.SubElement(position, "dZPosition", runtype="double", value=f"{z:.15f}")
    ET.SubElement(position, "dPFSOffset", runtype="double", value=f"{pfs_offset:.15f}")
    return position


def export_grid(rows, columns, delta_y, z, pfs_offset, preamble_elements):
    variant = ET.Element("variant", version="1.0")
    no_name = ET.SubElement(variant, "no_name", runtype="CLxListVariant")
    no_name.extend(preamble_elements)
    forward_row = True
    idx = 0
    for row_idx, row_y in enumerate(rows):
        columns_iter = enumerate(columns)
        if not forward_row:
            # enumerate isn't reversible, see https://github.com/cjrh/enumerate_reversible
            columns_iter = reversed(list(columns_iter))
        for col_idx, x in columns_iter:
            name = f"{row_idx:02d}.{col_idx:02d}"
            y = row_y + col_idx * delta_y
            no_name.append(position_element(idx, name, x, y, z, pfs_offset))
            idx += 1
        forward_row = not forward_row
    return ET.ElementTree(variant)


def autogrid_from_corners(input_xml, x_step, y_step, num_rows):
    if y_step is not None and num_rows is not None:
        raise ValueError("expecting either y_step or num_rows but not both")
    input_positions = []
    preamble_elements = []
    for e in input_xml.findall("./no_name/*"):
        if e.tag in ("bIncludeZ", "bPFSEnabled"):
            preamble_elements.append(e)
        else:
            input_positions.append(parse_position(e))
    if len(input_positions) != 4:
        raise ValueError("expecting four corners")
    upper_left = input_positions[0]
    upper_right = input_positions[1]
    lower_right = input_positions[2]
    lower_left = input_positions[3]
    # use z, pfs_offset from upper_left corner position
    z = upper_left["z"]
    pfs_offset = upper_left["pfs_offset"]
    if num_rows is not None:
        rows = np.linspace(upper_left["y"], lower_left["y"], num_rows)
    else:
        rows = np.arange(upper_left["y"], lower_left["y"], -y_step)
    columns = np.arange(upper_left["x"], upper_right["x"], -x_step)
    # delta_y is the increment in y every time you move to the next column
    slope = (upper_right["y"] - upper_left["y"]) / (upper_right["x"] - upper_left["x"])
    delta_y = -slope * x_step
    return export_grid(rows, columns, delta_y, z, pfs_offset, preamble_elements)


@click.group()
def cli():
    pass


@cli.command()
@click.argument("input", type=click.File("r", encoding="utf-16"))
# @click.argument("output", type=click.File("w", encoding="utf-16"))
@click.argument("output", type=click.Path())
@click.option("--x", type=float)  # these are given in µm
@click.option("--y", type=float)
@click.option("--rows", type=int)
def from_corners(input, output, x, y, rows):
    if y is not None and rows is not None:
        click.error("Expecting exactly one of --y and --rows")
    input_xml = ET.parse(input, parser=ET.XMLParser(encoding="utf-16"))
    output_xml = autogrid_from_corners(input_xml, x, y, rows)
    # this should write a UTF-16LE with BOM under windows
    # SEE: https://peter.bloomfield.online/why-python-3-doesnt-write-the-unicode-bom/
    output_xml.write(output, encoding="utf-16")


if __name__ == "__main__":
    cli()
