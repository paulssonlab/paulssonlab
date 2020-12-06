import numpy as np
import re
import requests
from copy import deepcopy
import string
import pygsheets
from Bio.Seq import Seq
import Bio.Restriction as Restriction
from paulssonlab.cloning.sequence import smoosh_sequences
from paulssonlab.api.addgene import get_addgene
from paulssonlab.api.google import (
    list_drive,
    upload_drive,
    columns_with_validation,
    insert_sheet_rows,
)
from paulssonlab.api import read_sequence, regex_key
from paulssonlab.api.util import PROGRESS_BAR

ID_REGEX = r"([A-Za-z]*)\s*(\d+(?:\.\d+)?)[a-zA-Z]*"

MARKER_ABBREVIATIONS = {
    "ampicillin": "Ampicillin/Carbenicillin (50 µg/mL)",
    "carbenicillin": "Ampicillin/Carbenicillin (50 µg/mL)",
    "kanamycin": "Kanamycin (50 µg/mL)",
    "chloramphenicol": "Chloramphenicol (25 µg/mL)",
    "tetracycline": "Tetracycline (10 µg/mL)",
    "spectinomycin": "Spectinomycin (50 µg/mL)",
    "gentamicin": "Gentamicin (20 µg/ml)",
}


def get_next_empty_row(worksheet, skip_columns=0):
    last_idx, _ = _get_next_empty_row(worksheet, skip_columns=skip_columns)
    if last_idx is None:
        return 2
    else:
        # increment twice for:
        # - add one to convert from zero-indexing to one-indexing
        # - row 1 is header
        return last_idx + 2


def _get_next_empty_row(worksheet, skip_columns=0):
    df = worksheet.get_as_df(has_header=False, empty_value=None)
    df = df[1:]
    has_datavalidation = columns_with_validation(
        worksheet.client.sheet.service,
        worksheet.spreadsheet.id,
        worksheet.spreadsheet._sheet_list,
    )
    mask = has_datavalidation[worksheet.title]
    mask += [False] * (len(df.columns) - len(mask))
    nonempty = (
        ~df.iloc[:, skip_columns:]
        .iloc[:, ~np.array(mask)[skip_columns:]]
        .isnull()
        .all(axis=1)
    )
    last_idx = nonempty[nonempty].last_valid_index()
    # convert to Python int because
    # DataFrame.last_valid_index() returns np.int64, which is not JSON-serializable
    last_idx = int(last_idx)
    return last_idx, df


def get_next_collection_id(worksheet):
    last_idx, df = _get_next_empty_row(worksheet, skip_columns=1)
    if last_idx is None:
        # sheet is empty, initialize at prefix 1
        prefix = worksheet.spreadsheet.title.split("_")[0]
        return (prefix, 1), 2
    # ID for last non-empty row
    last_id = df.iloc[last_idx - 1, 0]
    prefix, number = re.match(ID_REGEX, str(last_id)).groups()
    number_parts = number.split(".")
    if len(number_parts) == 2 and number_parts[0] == "0":
        # increment decimal if whole-part of the number is 0
        prefix += "0."
        index = int(number_parts[1])
    else:
        index = int(number_parts[0])
    # increment twice for:
    # - add one to convert from zero-indexing to one-indexing
    # - row 1 is header
    row = last_idx + 2
    return (prefix, index + 1), row


def trim_unassigned_ids(worksheet):
    _, row = get_next_collection_id(worksheet)
    _trim_unassigned_ids(worksheet, row)


def _trim_unassigned_ids(worksheet, row):
    values = [[""] * (worksheet.rows - row)]
    worksheet.update_values(f"A{row}:", values, majordim="COLUMNS")


def rename_ids(sheet, old_prefix, new_prefix, skiprows=1, column=1):
    values = sheet.get_values(
        (skiprows + 1, column),
        (sheet.rows, column),
        value_render=pygsheets.ValueRenderOption.FORMULA,
    )
    new_values = [
        [rename_id(value, old_prefix, new_prefix) for value in row] for row in values
    ]
    sheet.update_values((skiprows + 1, column), new_values, parse=True)


def rename_id(s, old_prefix, new_prefix):
    match = re.match(f"^{re.escape(old_prefix)}(.*)", s)
    if match:
        return f"{new_prefix}{match.group(1)}"
    else:
        try:
            float(s)
            return f"{new_prefix}{s}"
        except:
            # if s isn't a number nor begins with old_prefix,
            # this is unexpected, so give up and don't change the ID
            return s


def import_addgene(
    urls,
    strain_sheet,
    plasmid_sheet,
    plasmid_maps_folder,
    parts=False,
    strain_overrides=None,
    plasmid_overrides=None,
    strain_defaults={"Marker*": "Unknown"},
    plasmid_defaults=None,
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
        strain_defaults=strain_defaults,
        plasmid_defaults=plasmid_defaults,
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
                plasmid_maps[filename] = {
                    "content": content,
                    "mimetype": "chemical/seq-na-genbank",
                }
    # add strains to spreadsheet
    insert_sheet_rows(strain_sheet, strain_row, strains)
    # add plasmids to spreadsheet
    insert_sheet_rows(plasmid_sheet, plasmid_row, plasmids)
    # copy sequences to plasmid maps folder
    if service is None:
        service = strain_sheet.client.drive.service
    upload_plasmid_maps(service, plasmid_maps, plasmid_maps_folder, overwrite=overwrite)
    # TODO: add file metadata/row numbers to returned data?
    # trim extra rows
    if trim:
        trim_unassigned_ids(strain_sheet)
        trim_unassigned_ids(plasmid_sheet)
    return dict(strains=strains, plasmids=plasmids, plasmid_maps=plasmid_maps)


def upload_plasmid_maps(service, plasmid_maps, plasmid_maps_folder, overwrite=True):
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
        file_id = upload_drive(
            service,
            content=plasmid_map["content"],
            name=filename,
            file_id=file_id,
            mimetype=plasmid_map["mimetype"],
            parent=plasmid_maps_folder,
        )
        plasmid_map["id"] = file_id
    return plasmid_maps


def _import_addgene_data(
    urls,
    strain_overrides=None,
    plasmid_overrides=None,
    strain_defaults=None,
    plasmid_defaults=None,
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
                strain_defaults=strain_defaults,
                plasmid_defaults=plasmid_defaults,
                callback=callback,
            )
        )
    return data


def _format_addgene_for_spreadsheet(
    data,
    strain_overrides=None,
    plasmid_overrides=None,
    strain_defaults=None,
    plasmid_defaults=None,
    callback=None,
):
    if data["item"] == "Kit":
        entries = []
        for well in data["wells"]:
            entry = _format_addgene_for_spreadsheet(
                well,
                strain_overrides=strain_overrides,
                plasmid_overrides=plasmid_overrides,
                strain_defaults=strain_defaults,
                plasmid_defaults=plasmid_defaults,
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
        if strain_defaults:
            strain = {**strain_defaults, **strain}
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
                res = requests.get(seq_url)
                plasmid_map = read_sequence(res.content.decode("utf8"))
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
            if plasmid_defaults:
                plasmid = {**plasmid, **plasmid_defaults}
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


def insert_parts(part_sheet, parts, first_row):
    parts = deepcopy(parts)
    parts = _prepare_parts(part_sheet, parts, first_row)
    # part_sheet.spreadsheet.default_parse must be True for formula to be parsed, but this is True by default
    return insert_sheet_rows(part_sheet, first_row, parts)


def _prepare_parts(part_sheet, parts, first_row):
    cols = part_sheet.get_row(1)
    overhang1_col = string.ascii_uppercase[cols.index("Upstream overhang*")]
    overhang2_col = string.ascii_uppercase[cols.index("Downstream overhang*")]
    seq_col = string.ascii_uppercase[cols.index("Sequence*")]
    for enzyme_name in ("BsaI", "BsmBI", "AarI", "BbsI"):
        enzyme = getattr(Restriction, enzyme_name)
        site1 = enzyme.site
        site2 = str(Seq(site1).reverse_complement())
        for row, part in enumerate(parts, first_row):
            has_enzyme = f'=IF(AND(ISERROR(SEARCH("{site1}",${seq_col}{row})), ISERROR(SEARCH("{site2}",${seq_col}{row}))),"","yes")'
            part[enzyme_name] = has_enzyme
    for row, part in enumerate(parts, first_row):
        part_type = f"=ArrayFormula(INDEX('Part types'!A:A, MATCH(${overhang1_col}{row}&\" \"&${overhang2_col}{row},{{'Part types'!B:B&\" \"&'Part types'!C:C}},0)))"
        part["Type*"] = part_type
        correct_overhangs = f'=IF(AND(LEFT(${seq_col}{row},LEN(${overhang1_col}{row}))=${overhang1_col}{row},RIGHT(${seq_col}{row},LEN(${overhang2_col}{row}))=${overhang2_col}{row}),"","no")'
        part["Correct overhangs?"] = correct_overhangs
    return parts


def get_drive_seq(drive_service, root, name):
    return next(get_drive_seqs(drive_service, root, [name]))


def get_drive_seqs(drive_service, root, names):
    seq_files = list_drive(drive_service, root=root)
    for name in names:
        seq_file = regex_key(seq_files, f"{name}\\.", check_duplicates=True)["id"]
        seq = read_sequence(
            drive_service.files().get_media(fileId=seq_file).execute().decode("utf8")
        )
        yield seq


def add_overhangs(seq, overhangs):
    seq = smoosh_sequences(overhangs[0], seq)
    seq = smoosh_sequences(seq, overhangs[1])
    return seq


def add_flanks(seq, flanks):
    for flank_5prime, flank_3prime in flanks:
        seq = flank_5prime + seq + flank_3prime
    return seq
