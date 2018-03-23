import numpy as np
import click
from peewee import Model
from playhouse.apsw_ext import (
    APSWDatabase,
    TextField,
    IntegerField,
    DoubleField,
    BooleanField,
)
from json_util import JSONField
import pickle
import base64
import zarr
from datetime import datetime
from collections import OrderedDict, defaultdict
from contextlib import contextmanager
import os
import re
import gc
from IPython import embed
from util import tqdm_auto
from metadata import parse_nd2_file_metadata, parse_nikon_tiff_file_metadata

AGGREGATE_EXCLUDE = re.compile(
    r"^(Thumbs.*\.db|Pos7_TimelapseAnalysis\.mat|.*\.txt|\..*)$"
)
AGGREGATE_EXCLUDE_DIRS = re.compile(r"^ps$")

EXTENSION_ALIASES = {"tif": "tiff"}

METADATA_READERS = {
    "tiff": parse_nikon_tiff_file_metadata,
    "nd2": parse_nd2_file_metadata,
}

db = APSWDatabase(None)


class BaseModel(Model):
    class Meta:
        database = db


class Status(BaseModel):
    name = TextField(unique=True)


class File(BaseModel):
    path = TextField(unique=True)
    type = TextField()
    size = IntegerField()
    atime = DoubleField()
    mtime = DoubleField()
    ctime = DoubleField()
    aggregated = BooleanField()
    aggregated_path = TextField()
    metadata = JSONField()


def connect_db(db_file):
    db.init(db_file, pragmas=[("journal_mode", "wal")])
    db.connect()
    db.create_tables([Status, File])
    return db


@contextmanager
def _processing_step(name, reprocess=False):
    res = Status.select().where(Status.name == name)
    if not res or reprocess:
        Status.delete().where(Status.name == name)
        yield True
        Status.create(name=name)
    else:
        yield False


def _get_extension(path, aliases=EXTENSION_ALIASES):
    extension = path.split(".")[-1].lower()
    extension = aliases.get(extension, extension)
    return extension


def _scan_file(path):
    doc = {}
    doc["path"] = path
    stat = os.stat(path)
    doc["size"] = stat.st_size
    doc["atime"] = stat.st_atime
    doc["mtime"] = stat.st_mtime
    doc["ctime"] = stat.st_ctime
    # doc['atime'] = datetime.fromtimestamp(stat.st_atime)
    # doc['mtime'] = datetime.fromtimestamp(stat.st_mtime)
    # doc['ctime'] = datetime.fromtimestamp(stat.st_ctime)
    return doc


def _abbreviate_filename(filename, length, sep="..."):
    return filename[: length // 2 - len(sep)] + sep + filename[-length // 2 :]


def _scan_file_list(file_list, valid_extensions, extension_aliases=EXTENSION_ALIASES):
    pbar = tqdm_auto(file_list, desc="scanning for image files")
    num_files_to_process = 0
    skipped = []
    for path in pbar:
        pbar.set_postfix(to_process=num_files_to_process)
        path = path.strip()
        if not os.path.exists(path):
            skipped.append(path)
            continue
        root, file = os.path.split(path)
        extension = _get_extension(file)
        if extension in valid_extensions:
            yield path, extension, False
            num_files_to_process += 1
    print("skipped: {}".format(skipped))


# TODO: untested/unused
def _is_empty_excluding(root, dirs, exclude):
    for d in dirs:
        if not exclude.match(d):
            new_root = os.path.join(root, d)
            is_empty = _is_empty_excluding(new_root, os.listdir(new_root), exclude)
            if not is_empty:
                return False
    else:
        return True


def _should_aggregate_directory(root, dirs, extension, files_by_extension):
    if (
        len(files_by_extension[extension]) > 5
        and sum(
            len(files) for ext, files in files_by_extension.items() if ext != extension
        )
        / len(files_by_extension[extension])
        < 0.1
    ):
        return True
    else:
        return False


def _scan_directory(directory, valid_extensions, aggregate_extensions=None):
    pbar = tqdm_auto(os.walk(directory), desc="scanning for image files")
    num_files_to_process = 0
    skipped = []
    for root, dirs, files in pbar:
        pbar.set_postfix(to_process=num_files_to_process)
        # modifying dirs in place will affect which dirs are traversed
        # SEE: https://stackoverflow.com/questions/19859840/excluding-directories-in-os-walk
        dirs[:] = [
            d
            for d in dirs
            if not d.startswith(".")
            and not os.path.exists(os.path.join(root, d, ".zattrs"))
        ]
        files_by_extension = defaultdict(list)
        for file in files:
            extension = _get_extension(file)
            if extension in valid_extensions:
                path = os.path.join(root, file)
                files_by_extension[extension].append(path)
        for extension, paths in files_by_extension.items():
            if extension in aggregate_extensions and _should_aggregate_directory(
                root, dirs, extension, files_by_extension
            ):
                yield root, extension, paths
                num_files_to_process += 1
            else:
                for path in paths:
                    yield path, extension, False
                    num_files_to_process += 1
    print("skipped: {}".format(skipped))  # TODO: not used


@click.command()
@click.argument("in_path", type=click.Path(exists=True))
@click.argument(
    "out_path", type=click.Path(writable=True, dir_okay=False), required=False
)
@click.option("--rescan", is_flag=True, default=False)
@click.option(
    "--file-list", type=click.File("r")
)  # TODO: make this a flag, have file list as in_path instead?
@click.option("--skip", multiple=True, default=False)
@click.option("--aggregate", multiple=True, default=["tiff"])
@click.option("--metadata/--no-metadata", default=True)
@click.option("--duplicates/--no-duplicates", default=False)
def inventory(
    in_path, out_path, rescan, file_list, skip, aggregate, metadata, duplicates
):
    # use one unified files table
    # add stat information on initial scan
    # add metadata in decreasing order of size
    # duplicates: store hash of small crop
    if file_list is not None and rescan:
        raise click.fail("cannot specify both --file-list and --rescan")
    db = connect_db(out_path)
    try:
        ### SCAN
        valid_extensions = set(METADATA_READERS.keys())
        for extension in skip:
            extension = EXTENSION_ALIASES.get(extension.lower(), extension.lower())
            valid_extensions.remove(extension)
        with _processing_step("scan", reprocess=rescan) as process:
            if process:
                if file_list:
                    files_to_scan = _scan_file_list(file_list, valid_extensions)
                else:
                    files_to_scan = _scan_directory(
                        in_path, valid_extensions, aggregate_extensions=aggregate
                    )
                for path, extension, aggregated_files in files_to_scan:
                    if aggregated_files:
                        # TODO: pick oldest file as exemplar??
                        scanned_files = [_scan_file(f) for f in aggregated_files]
                        scanned_files = sorted(scanned_files, key=lambda x: x["mtime"])
                        total_size = sum(d["size"] for d in scanned_files)
                        doc = scanned_files[0]
                        doc["size"] = total_size
                        doc["type"] = extension
                        doc["metadata"] = ""
                        doc["aggregated"] = True
                        doc["aggregated_path"] = path
                        del scanned_files  # TODO: probably not necessary, GC will get it
                    else:
                        doc = _scan_file(path)
                        doc["type"] = extension
                        doc["metadata"] = ""
                        doc["aggregated"] = False
                        doc["aggregated_path"] = ""
                    File.create(**doc)
        ### METADATA
        if metadata:
            skipped = []
            files_to_process = (
                File.select().where(File.metadata == "").order_by(File.size.desc())
            )  # sorted(files_table.search(~Query().metadata.exists()), key=lambda x: x['size'], reverse=True)
            pbar = tqdm_auto(files_to_process)
            for file_row in pbar:
                path = file_row.path
                file = os.path.split(path)[-1]
                extension = _get_extension(file)
                pbar.set_postfix(
                    file=_abbreviate_filename(file, 20), skipped=len(skipped)
                )
                try:
                    file_row.metadata = METADATA_READERS[extension](path)
                except KeyboardInterrupt:
                    print("skipped: {}".format(skipped))
                    raise
                except:
                    skipped.append(path)
                else:
                    file_row.save()
                #     files_table.write_back([doc])
                #     gc.collect()
            print("skipped: {}".format(skipped))
    finally:
        db.close()  # flush


@click.command()
@click.argument("file", type=click.Path(exists=True, dir_okay=False))
@click.option("--force/--no-force", default=False)
def inspect_metadata(file, force):
    # convenience imports
    from metadata import _nikon_tiff_label, _nikon_tiff_field, parse_nikon_tiff_metadata
    import PIL.Image
    import nd2reader

    extension = _get_extension(os.path.split(file)[-1])
    if force:
        try:
            md = METADATA_READERS[extension](file)
        except:
            print("could not load metadata")
    else:
        md = METADATA_READERS[extension](file)
    embed()


commands = [inventory, inspect_metadata]
