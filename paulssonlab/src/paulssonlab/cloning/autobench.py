#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd
import sys
import os
import toml
from collections import defaultdict
from prompt_toolkit import prompt
import click
from prompt_toolkit import prompt
from prompt_toolkit.completion import Completer, Completion
from prompt_toolkit.validation import Validator, ValidationError
from IPython import embed
from functools import partial, wraps
import requests
import operator
from itertools import chain

BENCHLING_URI = "https://benchling.com/api/v2"
BENCHLING1_URI = "https://benchling.com/api/v1"
FORWARD_STRAND = 1
REVERSE_STRAND = -1
VARIANT_SEP = "+"
PART_SEP = "-"

default_config_file = os.path.join(os.path.dirname(__file__), "config.toml")

# def list_benchling_projects(benchling):
#     pass


@click.group()
@click.option("--config", type=click.Path(), default=default_config_file)
@click.option("--benchling-key", default=None)
@click.pass_context
def cli(ctx, config, benchling_key):
    if not config:
        config = default_config_file
    cfg = toml.load(config)
    if benchling_key:
        cfg["benchling_key"] = benchling_key
    elif not cfg["benchling_key"]:
        print("Need Benchling API key.")
        sys.exit(1)
    # cfg['benchling'] = BenchlingAPI(cfg['benchling_key'])
    ctx.obj = cfg


def apply_first(func_to_apply):
    def decorator(func):
        @wraps(func)
        def f(first, *args, **kwargs):
            return func(func_to_apply(first), *args, **kwargs)

        return f

    return decorator


def apply_prefix(prefix):
    return apply_first(lambda x: prefix + x)


def paginator(req_func):
    @wraps(req_func)
    def f(*args, **kwargs):
        res = req_func(*args, **kwargs).json()
        non_token_keys = [k for k in res if k != "nextToken"]
        if len(non_token_keys) != 1:
            raise Exception(
                "got unexpected number of keys in paginator: {}".format(res.keys())
            )
        key = non_token_keys[0]
        for item in res[key]:
            yield item
        while "nextToken" in res and res["nextToken"]:
            kwargs["params"] = {
                "nextToken": res["nextToken"],
                **(kwargs["params"] if "params" in kwargs else {}),
            }
            res = req_func(*args, **kwargs).json()
            for item in res[key]:
                yield item

    return f


PARTS = {
    "Bujard_RBS": "GAATTCATTAAAGAGGAGAAAGGTCAT",
    "proD": "CACAGCTAACACCACGTCGTCCCTATCTGCTGCCCTAGGTCTATGAGTGGTTGCTGGATAACTTTACGGGCATGCATAAGGCTCGTATAATATATTCAGGGAGACCACAACGGTTTCCCTCTACAAATAATTTTGTTTAACTT",
    "ITS4": "TCCTCCGCTTATTGATATGC",
    "mScarlet_I": "AACGGTCACGAGTTCGAGATCGAAGGCGAAGGCGAGGGCCGTCCGTATGAAGGCACCCAGACCGCCAAACTGAAAGTGACTAAAGGCGGCCCGCTGCCTTTTTCCTGGGACATCCTGAGCCCGCAATTTATGTACGGTTCTAGGGCGTTCATCAAACACCCAGCGGATATCCCGGACTATTATAAGCAGTCTTTTCCGGAAGGTTTCAAGTGGGAACGCGTAATGAATTTTGAAGATGGTGGTGCCGTGACCGTCACTCAGGACACCTCCCTGGAGGATGGCACCCTGATCTATAAAGTTAAACTGCGTGGTACTAATTTTCCACCTGATGGCCCGGTGATGCAGAAAAAGACGATGGGTTGGGAGGCGTCTACCGAACGCTTGTATCCGGAAGATGGTGTGCTGAAAGGCGACATTAAAATGGCCCTGCGCCTGAAAGATGGCGGCCGCTATCTGGCTGACTTCAAAACCACGTACAAAGCCAAGAAACCTGTGCAGATGCCTGGCGCGTACAATGTGGACCGCAAACTGGACATCACCTCTCATAATGAAGATTATACGGTGGTAGAGCAATATGAGCGCTCCGAGGGTCGTCATTCTACCGGTGGCATGGATGAACTATACAAA",
    "mSCFP3": "AGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTCTCACTTGGGGTGTTCAATGCTTTGCAAGATACCCAGATCATATGAAACAGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTACATCTCAGACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTAGACACAACATTGAAGATGGAGGCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTAAGCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAA",
    "rcVenus": "AGTAAAGGCGAAGAATTGTTCACTGGCGTGGTACCGATCCTGGTAGAACTGGATGGCGACGTTAATGGTCACAAGTTCAGCGTTAGTGGAGAGGGTGAAGGTGATGCGACCTATGGCAAACTGACCCTGAAGCTGATCTGCACAACCGGCAAGCTGCCTGTTCCTTGGCCGACACTGGTTACAACGCTGGGCTATGGCGTACAATGTTTCGCACGGTACCCGGACCACATGAAGCAACATGACTTCTTCAAGAGCGCTATGCCTGAAGGCTATGTCCAAGAAAGGACTATCTTCTTCAAAGACGACGGCAATTACAAGACACGGGCCGAAGTCAAATTCGAAGGCGATACGCTGGTCAACAGAATCGAGCTGAAAGGCATCGACTTCAAGGAAGATGGCAACATCCTGGGCCATAAACTGGAATATAATTATAACAGTCATAATGTGTATATCACCGCTGACAAACAAAAGAATGGCATCAAGGCCAACTTCAAAATCAGACATAACATCGAAGATGGAGGTGTTCAACTGGCAGACCACTACCAACAAAATACTCCGATCGGCGATGGCCCGGTGCTGCTGCCGGATAACCATTATCTGAGTTATCAAAGTAAGCTGAGCAAGGATCCGAACGAAAAAAGAGATCATATGGTTCTGCTGGAATTCGTAACGGCCGCGGGCATCACGCATGGCATGGACGAGCTGTATAAA",
    "Cluzel_RBS": "TAGTAAGGAGTCTAGACC",
}
PARTS = {k: v.lower() for k, v in PARTS.items()}
# TODO: grab parts list from benchling (use registry entry to mark parts?)


def find_parts(bases):
    for name, part in PARTS.items():
        idx = bases.find(part)
        if idx != -1:
            head = bases[:idx]
            tail = bases[idx + len(part) :]
            return (
                find_parts(head)
                + [
                    name,
                ]
                + find_parts(tail)
            )
    if len(bases):
        return [(bases,)]
    else:
        return []


def _format_part(part, variant_sep=VARIANT_SEP):
    if isinstance(part, tuple):
        if len(part) == 1:  # raw bases
            return "({})".format(len(part[0]))
        elif len(part) == 2:  # part variant
            return variant_sep.join([part[0]] + list(part[1]))
        else:
            raise ValueError("unexpected length of tuple part")
    else:
        return part


def format_parts(parts, sep="-"):
    return sep.join(map(_format_part, parts))


def summarize_parts(bases, sep="-"):
    return format_parts(find_parts(bases), sep=sep)


# find overlapping/similarly named primers for point mutations
# find f/r primers for amplification
def find_primer_pairs(seq):
    bases = seq["bases"]
    primers = seq["primers"]
    primers = sorted(primers, key=operator.itemgetter("start"))
    sections = []
    for idx, primer_f in enumerate(primers):
        if primer_f["strand"] == FORWARD_STRAND:
            primer_r = None
            for primer2 in chain(primers[idx + 1 :], primers[:idx]):
                if primer2["strand"] == REVERSE_STRAND:
                    primer_r = primer2
                    break
            if primer_r is None:
                print("could not find matching backwards primer")
            start = primer_f["start"]
            stop = primer_r["end"]
            if stop < start:
                section_bases = bases[start:] + bases[:stop]
            else:
                section_bases = bases[start:stop]
            sections.append(
                {
                    "start": start,
                    "stop": stop,
                    "length": len(section_bases),
                    "bases": section_bases,
                    "primer_f": primer_f,
                    "primer_r": primer_r,
                }
            )
    # _, backbone_idx = max(enumerate([s['length'] for s in sections]), key=operator.itemgetter(1))
    return sections


def _iter_folders_with_prefix(get, get1, prefixes):
    paginate = paginator(get)
    for folder in paginate("/folders"):
        if any(folder["name"].startswith(prefix) for prefix in prefixes):
            # TODO: why aren't we paginating?
            # for dna in paginate1('/sequences', params={'folder': folder['id']}):
            #     print('SEQ: ', dna['name'])
            res = get1("/sequences/", params={"folder": folder["id"]}).json()
            for seq in res["sequences"]:
                yield folder["name"], seq


def parse_variants(variant, variant_sep=VARIANT_SEP, part_sep=PART_SEP):
    parts = variant.split(part_sep)
    split_parts = [p.split(variant_sep) for p in parts]
    return [(p[0], p[1:]) for p in split_parts]


def strip_variants(variant, variant_sep=VARIANT_SEP, part_sep=PART_SEP):
    return part_sep.join(
        [
            v[0]
            for v in parse_variants(variant, variant_sep=variant_sep, part_sep=part_sep)
        ]
    )


def _iter_constructs(get, constructs):
    for construct_variant in constructs:
        construct = strip_variants(construct_variant)
        res = get("/dna-sequences", params={"name": construct})
        res.raise_for_status()
        seqs = res.json()["dnaSequences"]
        if len(seqs) == 0:
            print("did not find construct {}".format(construct))
            continue
        elif len(seqs) > 1:
            raise Exception("got more than one sequence with name {}".format(construct))
        yield construct_variant, seqs[0]


@cli.command()
@click.pass_obj
@click.argument("constructs", nargs=-1)
@click.option("--folder/--construct", default=True)
def show_primer_pairs(cfg, constructs, folder):
    s = requests.Session()
    s.auth = (cfg["benchling_key"], "")
    get = apply_prefix(BENCHLING_URI)(s.get)
    get1 = apply_prefix(BENCHLING1_URI)(s.get)
    paginate = paginator(get)
    paginate1 = paginator(get1)
    if folder:
        seqs = _iter_folders_with_prefix(get, get1, constructs)
    else:
        seqs = _iter_constructs(get, constructs)
    last_folder = None
    for name, seq in seqs:
        if folder:
            if name != last_folder:
                print("FOLDER: ", name)
                last_folder = name
        if not seq["primers"]:
            continue
        print("SEQ: ", seq["name"])
        print()
        sections = find_primer_pairs(seq)
        for section in sections:
            print("FWD: {}".format(section["primer_f"]["name"]))
            print(summarize_parts(section["bases"].lower()))
            print("REV: {}".format(section["primer_r"]["name"]))
            print("...")
        print()
        print("---")
        print()


@cli.command()
@click.pass_obj
@click.argument("constructs", nargs=-1)
@click.option("--folder/--construct", default=False)
@click.option("--start", default=1)
@click.option("--abbreviate-primers/--no-abbreviate-primers", default=True)
@click.option("--abbreviate-templates/--no-abbreviate-templates", default=True)
def plan_cloning(
    cfg, constructs, folder, start, abbreviate_primers, abbreviate_templates
):
    s = requests.Session()
    s.auth = (cfg["benchling_key"], "")
    get = apply_prefix(BENCHLING_URI)(s.get)
    get1 = apply_prefix(BENCHLING1_URI)(s.get)
    paginate = paginator(get)
    paginate1 = paginator(get1)
    if folder:
        seqs = _iter_folders_with_prefix(get, get1, constructs)
    else:
        seqs = _iter_constructs(get, constructs)
    transformations = []
    for construct, seq in seqs:
        if not seq["primers"]:
            continue
        variants = parse_variants(construct)
        sections = find_primer_pairs(seq)
        pcrs = []
        for section in sections:
            parts = find_parts(section["bases"].lower())
            for part, variant in variants:
                if variant:
                    parts = [(p, frozenset(variant)) if p == part else p for p in parts]
            pcrs.append(
                (section["primer_f"]["name"], section["primer_r"]["name"], tuple(parts))
            )
        transformations.append({"name": construct, "pcrs": pcrs})
    pcr_number = {}
    entries = []
    number = start
    for transformation in transformations:
        new_pcrs = [pcr for pcr in transformation["pcrs"] if pcr not in pcr_number]
        if new_pcrs:
            for pcr in new_pcrs:
                pcr_number[pcr] = number
                primer_f, primer_r, template = pcr
                if abbreviate_templates:
                    template = [
                        t for t in template if not isinstance(t, tuple) or len(t) != 1
                    ]
                template = format_parts(template)
                if abbreviate_primers:
                    primer_f = primer_f.split("_")[0]
                    primer_r = primer_r.split("_")[0]
                entry = {
                    "number": number,
                    "primer_f": primer_f,
                    "primer_r": primer_r,
                    "template": template,
                    "uses": sum(1 for t in transformations if pcr in t["pcrs"]),
                }
                number += 1
                entries.append(entry)
        else:
            entry = {"number": number}
            number += 1
            entries.append(entry)
        entry["name"] = transformation["name"]
        entry["template plasmid"] = ""
        entry["ingredients"] = "+".join(
            [
                "THIS" if pcr_number[pcr] == entry["number"] else str(pcr_number[pcr])
                for pcr in transformation["pcrs"]
            ]
        )
    spreadsheet = pd.DataFrame(entries).T.reindex(
        [
            "number",
            "name",
            "ingredients",
            "uses",
            "template",
            "template_plasmid",
            "primer_f",
            "primer_r",
        ]
    )
    print(spreadsheet.to_csv(None, sep="\t", header=False, index=False))


if __name__ == "__main__":
    cli()
