def format_well_name(plate, row, column):
    return f"{plate if plate and plate != 1 else ''}{row}{column}"
