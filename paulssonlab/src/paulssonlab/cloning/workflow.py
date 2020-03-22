import numpy as np
import re
from paulssonlab.api.google import (
    get_drive_by_name,
    filter_drive,
    columns_with_validation,
)


def get_strain_collection_sheets(service, collection_prefix):
    collection_folder = get_drive_by_name(
        service, f"{collection_prefix}_Collection", folder=True
    )
    files = (
        service.files()
        .list(q=f"'{collection_folder}' in parents")
        .execute()
        .get("files", [])
    )
    keys = {
        "strains": (f"{collection_prefix}_strains", False),
        "oligos": (f"o{collection_prefix}_oligos", False),
        "plasmids": (f"p{collection_prefix}_plasmids", False),
        "parts": (f"{collection_prefix}_parts", False),
        "plasmid_maps": (f"Plasmid_Maps", True),
    }
    collection = filter_drive(files, keys)
    return {"root": collection_folder, **collection}


def get_next_collection_id(worksheet):
    # ids = worksheet.get_values('A', 'A')
    # last_id = ids[-1][0]
    df = worksheet.get_as_df(has_header=False)
    df = df[1:]
    has_datavalidation = columns_with_validation(
        worksheet.client.sheet.service,
        worksheet.spreadsheet.id,
        worksheet.spreadsheet._sheet_list,
    )
    mask = has_datavalidation[worksheet.title]
    mask += [False] * (len(df.columns) - len(mask))
    nonempty = ~df.iloc[:, ~np.array(mask)].iloc[:, 1:].isnull().any(axis=1)
    last_idx = nonempty[nonempty].last_valid_index()
    if last_idx is None:
        last_idx = len(nonempty)
    last_idx -= 1
    last_id = df.iloc[last_idx, 0]
    prefix, index, _ = re.match(r"([A-Za-z]+)(\d+)(\.\d+\w+?)?", last_id).groups()
    row_num = last_idx + 3  # increment, header, and one-indexing
    return (prefix, int(index) + 1), row_num
