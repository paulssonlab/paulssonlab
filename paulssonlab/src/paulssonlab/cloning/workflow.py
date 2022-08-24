import numpy as np
import re
import requests
from copy import deepcopy
import string
from datetime import datetime
from collections import ChainMap
import textwrap
import pygsheets
from Bio.Seq import Seq
import Bio.Restriction as Restriction
from paulssonlab.cloning.sequence import (
    get_seq,
    smoosh_sequences,
    find_homologous_ends,
    find_aligned_subsequence,
    DsSeqRecord,
    reverse_complement,
    assemble,
)
from paulssonlab.cloning.commands.parser import expr_list_parser
from paulssonlab.cloning.enzyme import re_digest
from paulssonlab.api.addgene import get_addgene
from paulssonlab.api.google import (
    list_drive,
    upload_drive,
    columns_with_validation,
    insert_sheet_rows,
)
from paulssonlab.api.util import PROGRESS_BAR, regex_key
from paulssonlab.cloning.io import bytes_to_value, filename_to_mimetype

DEGENERATE_BASES = "RYMKSWHBVDN".lower()
DEGENERATE_BASES_REGEX = re.compile(f"[{DEGENERATE_BASES}]", re.IGNORECASE)

ID_REGEX = r"\s*([A-Za-z]*)\s*(\d+)(?:[.a-zA-Z]|\s+|$)\S*\s*$"

# these are not currently used, but can be used for formatted outputs
# these numbers are taken from Addgene: https://www.addgene.org/mol-bio-reference/
# (see also https://blog.addgene.org/plasmids-101-everything-you-need-to-know-about-antibiotic-resistance-genes)
MARKER_DEFAULT_CONCENTRATIONS = {
    "amp": "Ampicillin/Carbenicillin (100 µg/mL)",
    "kan": "Kanamycin (50 µg/mL)",
    "chlor": "Chloramphenicol (25 µg/mL)",
    "tet": "Tetracycline (10 µg/mL)",
    "spec": "Spectinomycin (50 µg/mL)",
}

MARKER_ABBREVIATIONS = {
    "ampicillin": "amp",
    "carbenicillin": "amp",
    "kanamycin": "kan",
    "chloramphenicol": "chlor",
    "tetracycline": "tet",
    "spectinomycin": "spec",
}


def parse_id(s):
    match = re.match(ID_REGEX, s)
    if match is None:
        raise ValueError(f"could not parse ID: '{s}'")
    prefix = match.group(1)
    index = int(match.group(2))
    return prefix, index


def format_abx_marker(s):
    marker = s.lower()
    # we remove commas, because we use commas to indicate multiple markers
    marker = marker.replace(",", "")
    # we match ug/ instead of ug in case an abx name contains “ug”
    marker = re.sub(r"ug\s*/\s*", "μg/", marker, flags=re.IGNORECASE)
    # capitalization of mL
    marker = marker.replace("ml", "mL")
    for name, abbrev in MARKER_ABBREVIATIONS.items():
        marker = re.sub(name, abbrev, marker, flags=re.IGNORECASE)
    return marker


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
    strain_defaults={"Marker": "Unknown"},
    plasmid_defaults=None,
    callback=None,
    trim=True,
    service=None,
    overwrite=True,
    progress_bar=PROGRESS_BAR,
):
    raise NotImplementedError  # TODO: this needs to be updated to new registry API
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
            source = entry[entry_type]["Source"]
            if source:
                components = source.split()
                if re.match(r"^https?://|www\.", components[0]):
                    source = f'=HYPERLINK("{components[0]}","{source}")'
                    entry[entry_type]["Source"] = source
        strain_id = f"{strain_prefix}{strain_number}"
        plasmid_id = f"{plasmid_prefix}{plasmid_number}"
        strain = entry.get("strain")
        if strain:
            strains.append(strain)
            strain["ID"] = strain_id
            strain_number += 1
            if "plasmid" in entry:
                strain["Plasmids"] = plasmid_id
        plasmid = entry.get("plasmid")
        if plasmid:
            plasmids.append(plasmid)
            plasmid["ID"] = plasmid_id
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
    # if trim:
    #     trim_unassigned_ids(strain_sheet)
    #     trim_unassigned_ids(plasmid_sheet)
    return dict(strains=strains, plasmids=plasmids, plasmid_maps=plasmid_maps)


def _import_addgene_data(
    urls,
    strain_overrides=None,
    plasmid_overrides=None,
    strain_defaults=None,
    plasmid_defaults=None,
    callback=None,
    progress_bar=PROGRESS_BAR,
):
    if isinstance(urls, (str, int)):
        urls = [urls]
    data = []
    if progress_bar is not None and len(urls) >= 2:
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
                entry["strain"]["Source"] += kit_source
            if "plasmid" in entry:
                entry["plasmid"]["Source"] += kit_source
            entries.append(entry)
        return entries
    elif data["item"] in ("Plasmid", "Bacterial Strain"):
        # grab copy number from plasmid map, fall back on copy number
        # put host strain in strain
        # store vector backbone in (???)
        background = data["growth strain(s)"]
        marker = format_abx_marker(data["bacterial resistance(s)"])
        source = data["url"]
        reference = data["how_to_cite"].get("references")
        strain = {
            "Names": data["name"],
            "Species": "E. coli",
            # "Genotype": "",
            "Background": background,
            "Parent": background,
            "Marker": marker,
            # "DNA Transformed": "",
            "Source": source,
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
                plasmid_map = bytes_to_value(res.content, filename_to_mimetype(seq_url))
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
                "Size": size,
                "Origin": origin,
                "Marker": marker,
                "Other Notes": other_notes,
                "Source": source,
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


def part_entry_to_seq(entry):
    seq = entry["Sequence"]
    upstream_overhang_seq = entry["Upstream overhang"]
    downstream_overhang_seq = entry["Downstream overhang"]
    if upstream_overhang_seq.lower() != seq[: len(upstream_overhang_seq)].lower():
        raise ValueError("upstream overhang does not match part sequence")
    if (
        downstream_overhang_seq.lower()
        != seq[len(seq) - len(downstream_overhang_seq) :].lower()
    ):
        raise ValueError("downstream overhang does not match part sequence")
    # TODO: this assumes a particular overhang convention for golden gate/type IIS restriction enzymes
    upstream_overhang = len(upstream_overhang_seq)
    downstream_overhang = -len(downstream_overhang_seq)
    return DsSeqRecord(
        Seq(seq),
        upstream_overhang=upstream_overhang,
        downstream_overhang=downstream_overhang,
    )


def re_digest_part(seq, enzyme):
    # TODO: this doesn't handle the case where there is a cut site in the backbone
    # not sure there's a general way to pick out the intended part in that case
    frags = re_digest(seq, enzyme)
    # TODO: this relies on inward_cut=None for blunt end cuts
    # we should probably make this more robust by explictly checking for upstream/downstream_overhang==0
    non_backbone = [
        f
        for f in frags
        if not (
            f.upstream_inward_cut is not True and f.downstream_inward_cut is not True
        )
    ]
    if len(non_backbone) == 0:
        raise ValueError(
            f"no fragments with only outward cuts were generated by {enzyme} digest"
        )
    product = assemble(non_backbone, method="goldengate")
    return product


def is_bases(s):
    return re.match(r"^[atcgrymkswbdhvn]+$", s, re.IGNORECASE) is not None


def concatenate_flanks(*flanks, lower=True):
    upstream, downstream = flanks[0]
    if len(flanks) >= 1:
        for upstream_flank, downstream_flank in flanks[1:]:
            upstream = upstream_flank + upstream
            downstream = downstream + downstream_flank
    if lower:
        upstream = upstream.lower()
        downstream = downstream.lower()
    return upstream, downstream


def smoosh_flanks(*flanks, lower=True):
    upstream, downstream = flanks[0]
    if lower:
        upstream = upstream.lower()
        downstream = downstream.lower()
    if len(flanks) >= 1:
        for upstream_flank, downstream_flank in flanks[1:]:
            if lower:
                upstream_flank = upstream_flank.lower()
                downstream_flank = downstream_flank.lower()
            upstream = smoosh_sequences(upstream_flank, upstream)
            downstream = smoosh_sequences(downstream, downstream_flank)
    return upstream, downstream


def smoosh_and_trim_flanks(seq, flanks, lower=True):
    if lower:
        seq = seq.lower()
    upstream_overlap = find_homologous_ends(flanks[0], seq)
    downstream_overlap = find_homologous_ends(seq, flanks[1])
    return (
        flanks[0][: len(flanks[0]) - upstream_overlap],
        flanks[1][downstream_overlap:],
    )


# TODO: update the above to normalize similarly?
def smoosh_and_normalize_sequences(*seqs):
    seqs = [normalize_seq(seq) for seq in seqs]
    return smoosh_sequences(*seqs)


def find_coding_sequence(seq, prefix="atg", suffix=["taa", "tga", "tag"]):
    seq_str = str(get_seq(seq)).lower()
    start = seq_str.find(prefix)
    stop = -1
    for suffix_ in suffix:
        new_stop = find_aligned_subsequence(seq_str[start:], suffix_, last=True)
        if new_stop is not None:
            new_stop += len(suffix_)
            if new_stop > stop:
                stop = new_stop
    if stop == -1:
        raise ValueError(f"could not find aligned suffix {suffix}")
    # stop is indexed from start, see call to find_aligned_subsequence above
    stop += start
    return start, stop


def get_source_plasmid(registry, usage):
    cmds = expr_list_parser.parse(usage)
    for cmd in cmds:
        if cmd["_type"] == "digest" and cmd["input"]["_type"] == "name":
            name = cmd["input"]["name"]
            prefix, types = registry.get_prefix_and_types(name)
            if "plasmids" in types:
                return name
    return None


def overhangs_for(x):
    return (
        normalize_seq(x["Upstream overhang"]),
        normalize_seq(x["Downstream overhang"]),
    )


def part_types_map(part_types):
    d = {}
    for name, row in part_types.items():
        d[overhangs_for(row)] = row["Type"]
    return d


def normalize_seq(seq):
    return str(get_seq(seq)).lower()


def normalize_seq_upper(seq):
    return str(get_seq(seq)).upper()


def date():
    return datetime.now().strftime("%-m/%-d/%Y")


def cluster_sequences_by_prefix(d, start_index=0, num_subseqs=None, max_length=None):
    seqs = list(d.values())
    if num_subseqs is None:
        num_subseqs = min(len(seq) for seq in seqs)
    if max_length is None:
        max_length = max(len(subseq) for seq in seqs for subseq in seq[:num_subseqs])
    mismatch = False
    for idx in range(max_length):
        for subseq_idx in range(num_subseqs):
            base0 = seqs[0][subseq_idx][idx]
            if not all(seq[subseq_idx][idx] == base0 for seq in seqs[1:]):
                mismatch = True
                break
        if mismatch:
            break
    if not mismatch:
        # idx refers to the first mismatched base, so if no mismatches found, set to max_length
        idx += 1
    common_key = tuple(
        seqs[0][subseq_idx][start_index:idx] for subseq_idx in range(num_subseqs)
    )
    if idx > start_index:
        if idx == max_length:
            values = tuple(d.keys())
        else:
            values = cluster_sequences_by_prefix(
                d, start_index=idx, num_subseqs=num_subseqs, max_length=max_length
            )
        res = {common_key: values, "_size": len(d)}
        return res
    else:
        clusters = {}
        for name, seq in d.items():
            key = tuple(
                seq[subseq_idx][start_index : idx + 1]
                for subseq_idx in range(num_subseqs)
            )
            clusters.setdefault(key, {})
            clusters[key][name] = seq
        res = ChainMap(
            {"_size": len(d)},
            *[
                cluster_sequences_by_prefix(
                    cluster,
                    start_index=idx,
                    num_subseqs=num_subseqs,
                    max_length=max_length,
                )
                for cluster in clusters.values()
            ],
        )
        return res


def print_sequence_clusters(
    clusters, length=0, indent_level=0, extra_indent=4, wrap_width=100
):
    indent_str = " " * indent_level
    for key, cluster in clusters.items():
        if key == "_size":
            continue
        segment = "/".join(key)
        if isinstance(cluster, tuple):
            num_seqs = len(cluster)
        else:
            num_seqs = cluster["_size"]
        new_length = length + len(key[0])
        print(f"{indent_str}{segment} ({new_length}nt x {num_seqs}):")
        if isinstance(cluster, tuple):
            print(
                textwrap.fill(
                    ", ".join(cluster),
                    width=wrap_width,
                    initial_indent=" " * (indent_level + extra_indent),
                    subsequent_indent=" " * (indent_level + extra_indent),
                    break_long_words=True,
                )
            )
            print()
        else:
            print_sequence_clusters(
                cluster, length=new_length, indent_level=indent_level + 2
            )
