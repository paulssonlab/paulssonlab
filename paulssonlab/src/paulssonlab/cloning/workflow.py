import numpy as np
import re
from paulssonlab.api.addgene import get_addgene
from paulssonlab.api.google import (
    get_drive_by_name,
    filter_drive,
    list_drive,
    upload_drive,
    columns_with_validation,
)
from paulssonlab.api import get_genbank
from paulssonlab.api.util import PROGRESS_BAR


MARKER_ABBREVIATIONS = {
    "ampicillin": "Ampicillin/Carbenicillin (50 µg/mL)",
    "carbenicillin": "Ampicillin/Carbenicillin (50 µg/mL)",
    "kanamycin": "Kanamycin (50 µg/mL)",
    "chloramphenicol": "Chloramphenicol (25 µg/mL)",
    "tetracycline": "Tetracycline (10 µg/mL)",
    "spectinomycin": "Spectinomycin (50 µg/mL)",
    "gentamicin": "Gentamicin (20 µg/ml)",
}

DEFAULT_VALUES = {"Marker*": "Unknown"}


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
    df = worksheet.get_as_df(has_header=False, empty_value=None)
    df = df[1:]
    has_datavalidation = columns_with_validation(
        worksheet.client.sheet.service,
        worksheet.spreadsheet.id,
        worksheet.spreadsheet._sheet_list,
    )
    mask = has_datavalidation[worksheet.title]
    mask += [False] * (len(df.columns) - len(mask))
    nonempty = ~df.iloc[:, 1:].iloc[:, ~np.array(mask)[1:]].isnull().all(axis=1)
    last_idx = nonempty[nonempty].last_valid_index()
    if last_idx is None:
        # sheet is empty, initialize at prefix 1
        prefix = worksheet.spreadsheet.title.split("_")[0]
        return (prefix, 1), 2
    else:
        # convert to Python int because
        # DataFrame.last_valid_index() returns np.int64, which is not JSON-serializable
        last_idx = int(last_idx)
    last_idx -= 1
    last_id = df.iloc[last_idx, 0]
    prefix, index, _ = re.match(r"([A-Za-z]*)(\d+)(\.\d+\w+?)?", str(last_id)).groups()
    # increment three times for:
    # - one after last strain number
    # - one row taken up by header
    # - add one to convert from zero-indexing to one-indexing
    row = last_idx + 3
    return (prefix, int(index) + 1), row


def trim_unassigned_ids(worksheet):
    _, row = get_next_collection_id(worksheet)
    _trim_unassigned_ids(worksheet, row)


def _trim_unassigned_ids(worksheet, row):
    values = [[""] * (worksheet.rows - row)]
    worksheet.update_values(f"A{row}:", values, majordim="COLUMNS")


def _insert_rows(sheet, row, entries, default_values):
    columns = sheet.get_row(1)
    values = [
        [entry.get(col, default_values.get(col, "")) for col in columns]
        for entry in entries
    ]
    # we insert at row - 1 because this inserts below the given row number
    sheet.insert_rows(row - 1, number=len(values), values=values)


def import_addgene(
    urls,
    strain_sheet,
    plasmid_sheet,
    plasmid_maps_folder,
    parts=False,
    strain_overrides=None,
    plasmid_overrides=None,
    callback=None,
    trim=True,
    service=None,
    overwrite=True,
    progress_bar=PROGRESS_BAR,
):
    data = _import_addgene_data(
        urls,
        progress_bar=progress_bar,
        strain_overrides=strain_overrides,
        plasmid_overrides=plasmid_overrides,
        callback=callback,
    )
    # assign LIB/pLIB numbers, add pLIB genotype to strain
    (strain_prefix, strain_number), strain_row = get_next_collection_id(strain_sheet)
    (plasmid_prefix, plasmid_number), plasmid_row = get_next_collection_id(
        plasmid_sheet
    )
    strains = []
    plasmids = []
    plasmid_maps = {}
    for entry in data:
        for entry_type in ("strain", "plasmid"):
            if entry_type not in entry:
                continue
            source = entry[entry_type]["Source*"]
            if source:
                components = source.split()
                if re.match(r"^https?://|www\.", components[0]):
                    source = f'=HYPERLINK("{components[0]}","{source}")'
                    entry[entry_type]["Source*"] = source
        strain_id = f"{strain_prefix}{strain_number}"
        plasmid_id = f"{plasmid_prefix}{plasmid_number}"
        strain = entry.get("strain")
        if strain:
            strains.append(strain)
            strain["ID*"] = strain_id
            strain_number += 1
            if "plasmid" in entry:
                strain["Plasmid(s)*"] = plasmid_id
        plasmid = entry.get("plasmid")
        if plasmid:
            plasmids.append(plasmid)
            plasmid["ID*"] = plasmid_id
            plasmid_number += 1
            plasmid_map = entry.get("plasmid_map")
            if plasmid_map:
                filename = f"{plasmid_id}.gbk"
                content = plasmid_map.format("genbank")
                plasmid_maps[filename] = {"content": content}
    # add entries to spreadsheets
    _insert_rows(strain_sheet, strain_row, strains, DEFAULT_VALUES)
    # add plasmids to spreadsheet
    _insert_rows(plasmid_sheet, plasmid_row, plasmids, DEFAULT_VALUES)
    # copy sequences to plasmid maps folder
    if service is None:
        service = strain_sheet.client.drive.service
    files = list_drive(service, root=plasmid_maps_folder)
    for filename, plasmid_map in plasmid_maps.items():
        if filename in files:
            if not overwrite:
                raise ValueError(
                    f"enable overwrite or remove existing plasmid map: {filename}"
                )
            file_id = files[filename]["id"]
        else:
            file_id = None
        mimetype = "chemical/seq-na-genbank"
        file_id = upload_drive(
            service,
            content=plasmid_map["content"],
            name=filename,
            file_id=file_id,
            mimetype=mimetype,
            parent=plasmid_maps_folder,
        )
        plasmid_map["id"] = file_id
    # TODO: add file metadata/row numbers to returned data?
    # trim extra rows
    if trim:
        trim_unassigned_ids(strain_sheet)
        trim_unassigned_ids(plasmid_sheet)
    return dict(strains=strains, plasmids=plasmids, plasmid_maps=plasmid_maps)


def _import_addgene_data(
    urls,
    strain_overrides=None,
    plasmid_overrides=None,
    callback=None,
    progress_bar=PROGRESS_BAR,
):
    if isinstance(urls, str):
        urls = [urls]
    data = []
    if progress_bar is not None and len(urls) > 1:
        urls = progress_bar(urls)
    for url in urls:
        addgene = get_addgene(
            url, include_sequences=True, include_details=True, progress_bar=progress_bar
        )
        data.extend(
            _format_addgene_for_spreadsheet(
                addgene,
                strain_overrides=strain_overrides,
                plasmid_overrides=plasmid_overrides,
                callback=callback,
            )
        )
    return data


def _format_addgene_for_spreadsheet(
    data, strain_overrides=None, plasmid_overrides=None, callback=None
):
    if data["item"] == "Kit":
        entries = []
        for well in data["wells"]:
            entry = _format_addgene_for_spreadsheet(
                well,
                strain_overrides=strain_overrides,
                plasmid_overrides=plasmid_overrides,
                callback=callback,
            )
            if len(entry):  # skip entries for which callback returned False
                entry = entry[0]
            else:
                continue
            kit_source = f" (from kit {data['url']})"
            if "strain" in entry:
                entry["strain"]["Source*"] += kit_source
            if "plasmid" in entry:
                entry["plasmid"]["Source*"] += kit_source
            entries.append(entry)
        return entries
    elif data["item"] in ("Plasmid", "Bacterial Strain"):
        # grab copy number from plasmid map, fall back on copy number
        # put host strain in strain
        # store vector backbone in (???)
        background = data["growth strain(s)"]
        marker = data["bacterial resistance(s)"].replace("μg/ml", "μg/mL")
        marker = MARKER_ABBREVIATIONS.get(marker.lower(), marker)
        source = data["url"]
        reference = data["how_to_cite"].get("references")
        strain = {
            "Names": data["name"],
            "Species*": "E. coli",
            #'Genotype*': '',
            "Background*": background,
            "Parent*": background,
            "Marker*": marker,
            #'DNA Transformed': '',
            "Source*": source,
            "Reference": reference,
        }
        other_notes = []
        if data.get("well"):
            other_notes.append(f"Well: {data['well']}")
        if data.get("growth instructions"):
            other_notes.append(f"Growth instructions: {data['growth instructions']}")
        if data.get("growth temperature") and "37" not in data["growth temperature"]:
            other_notes.append(f"Growth temperature: {data['growth temperature']}")
        if data.get("depositor comments"):
            other_notes.append(f"Depositor comments: {data['depositor comments']}")
        other_notes = "\n".join(other_notes) or None
        strain["Other Notes"] = other_notes
        if strain_overrides:
            strain = {**strain, **strain_overrides}
        if data["item"] == "Plasmid":
            seq_url = None
            for seq_type in ("depositor_full", "addgene_full"):
                seq_urls = data["sequence_urls"].get(seq_type)
                if not seq_urls:
                    continue
                if len(seq_urls) > 1:
                    raise ValueError(
                        f"expecting one {seq_type.replace('_', ' ')} sequence"
                    )
                elif len(seq_urls) == 1:
                    seq_url = seq_urls[0]
                    break
            if seq_url is None:
                plasmid_map = None
                size = None
                origin = data.get("copy number") or "Unknown"
            else:
                plasmid_map = get_genbank(seq_url)
                size = len(plasmid_map.seq)
                ori_feature = next(
                    f for f in plasmid_map.features if f.type == "rep_origin"
                )
                origin = ori_feature.qualifiers["label"][0].strip()
                if origin == "ori":
                    origin = ori_feature.qualifiers["note"][0]
                else:
                    origin = re.sub(r" (?:ori|origin)$", "", origin)
            plasmid = {
                "Names": data["name"],
                "Size (bp)": size,
                "Origin*": origin,
                "Marker*": marker,
                "Other Notes": other_notes,
                "Source*": source,
                "Reference": reference,
            }
            if data.get("purpose"):
                plasmid["Description"] = data["purpose"]
            if plasmid_overrides:
                plasmid = {**plasmid, **plasmid_overrides}
            entry = dict(strain=strain, plasmid=plasmid, plasmid_map=plasmid_map)
        else:
            if data.get("purpose"):
                strain["Description"] = data["purpose"]
            entry = dict(strain=strain)
        if callback:
            entry = callback(entry, data)
        if entry is False:  # if callback returns False, skip entry
            return []
        else:
            return [entry]
    else:
        raise ValueError(f"unknown Addgene item type: {data['item']}")
