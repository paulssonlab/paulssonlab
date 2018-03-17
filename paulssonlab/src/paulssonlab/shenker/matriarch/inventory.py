import numpy as np
import click
from peewee import Model
from playhouse.apsw_ext import APSWDatabase, TextField, IntegerField, DoubleField
from json_util import JSONField
import pickle
import base64
import zarr
from datetime import datetime
from collections import OrderedDict
from contextlib import contextmanager
import os
import gc
from util import tqdm_auto
from metadata import parse_nd2_file_metadata, parse_nikon_tiff_file_metadata

METADATA_READERS = {
    "tif": parse_nikon_tiff_file_metadata,
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


def _get_extension(path):
    return path.split(".")[-1].lower()


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


@click.command()
@click.argument("in_path", type=click.Path(exists=True))
@click.argument(
    "out_path", type=click.Path(writable=True, dir_okay=False), required=False
)
@click.option("--rescan", is_flag=True, default=False)
@click.option("--file-list", type=click.File("r"))
@click.option("--skip-tiff", is_flag=True, default=False)
@click.option("--skip-nd2", is_flag=True, default=False)
@click.option("--metadata/--no-metadata", default=True)
@click.option("--duplicates/--no-duplicates", default=False)
def inventory(
    in_path, out_path, rescan, file_list, skip_tiff, skip_nd2, metadata, duplicates
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
        if skip_tiff:
            valid_extensions -= set(["tif", "tiff"])
        if skip_nd2:
            valid_extensions.remove("nd2")
        with _processing_step("scan", reprocess=rescan) as process:
            if process:
                skipped = []
                num_files_to_process = 0
                if file_list:
                    pbar = tqdm_auto(file_list, desc="scanning for image files")
                    for path in pbar:
                        path = path.strip()
                        if not os.path.exists(path):
                            skipped.append(path)
                            continue
                        root, file = os.path.split(path)
                        extension = _get_extension(file)
                        if extension in valid_extensions:
                            doc = _scan_file(path)
                            doc["type"] = extension
                            doc["metadata"] = ""
                            File.create(**doc)
                            num_files_to_process += 1
                else:
                    pbar = tqdm_auto(os.walk(in_path), desc="scanning for image files")
                    for root, dirs, files in pbar:
                        pbar.set_postfix(to_process=num_files_to_process)
                        # modifying dirs in place will affect which dirs are traversed
                        # SEE: https://stackoverflow.com/questions/19859840/excluding-directories-in-os-walk
                        dirs[:] = [
                            d
                            for d in dirs
                            if not os.path.exists(os.path.join(root, d, ".zattrs"))
                        ]
                        for file in files:
                            extension = _get_extension(file)
                            if extension in valid_extensions:
                                path = os.path.join(root, file)
                                if not os.path.exists(path):
                                    skipped.append(path)
                                    continue
                                doc = _scan_file(path)
                                doc["type"] = extension
                                doc["metadata"] = ""
                                File.create(**doc)
                                num_files_to_process += 1
                print("skipped: {}".format(skipped))
        ### METADATA
        with _processing_step("metadata") as process:
            if process:
                skipped = []
                files_to_process = (
                    File.select()
                    .where(File.metadata.is_null(False))
                    .order_by(File.size.desc())
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
