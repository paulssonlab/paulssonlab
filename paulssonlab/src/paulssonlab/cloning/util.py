import string
from itertools import count, product


def well_iterator(kind=96):
    if kind != 96:
        raise NotImplementedError
    for col, row in product(range(8), range(1, 13)):
        yield format_well_name(1, string.ascii_uppercase[col], row)


def format_well_name(plate_num, row, column):
    return f"{plate_num if plate_num and plate_num != 1 else ''}{row}{column}"
