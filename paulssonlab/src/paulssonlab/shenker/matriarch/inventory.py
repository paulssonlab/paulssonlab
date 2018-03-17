import numpy as np
import click
from tinydb import TinyDB, Query
from tinydb.middlewares import CachingMiddleware
from tinydb_serialization import SerializationMiddleware, Serializer
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

# FROM: https://github.com/msiemens/tinydb-serialization
# FROM: https://github.com/msiemens/tinydb-serialization/issues/6
class DateTimeSerializer(Serializer):
    OBJ_CLASS = datetime  # The class this serializer handles

    def encode(self, obj):
        return pickle.dumps(obj)  # obj.strftime('%Y-%m-%dT%H:%M:%S.%f')

    def decode(self, s):
        return pickle.loads(s)  # datetime.strptime(s, '%Y-%m-%dT%H:%M:%S.%f')


class ZarrSerializer(Serializer):
    OBJ_CLASS = zarr.Array

    def encode(self, obj):
        # return pickle.dumps(obj)
        return base64.encodebytes(pickle.dumps(obj)).decode()

    def decode(self, s):
        # return pickle.loads(s)
        return pickle.loads(base64.decodebytes(s.encode()))


class BytesSerializer(Serializer):
    OBJ_CLASS = bytes

    def encode(self, obj):
        return base64.encodebytes(obj).decode()

    def decode(self, s):
        return base64.decodebytes(s.encode())


def get_tinydb(db_file):
    serialization = SerializationMiddleware()
    # serialization.register_serializer(DateTimeSerializer(), 'datetime')
    serialization.register_serializer(ZarrSerializer(), "zarr")
    serialization.register_serializer(BytesSerializer(), "bytes")
    storage = CachingMiddleware(serialization)
    storage.WRITE_CACHE_SIZE = 50  # flush every 50 writes
    db = TinyDB(db_file, storage=storage)
    return db


def _get_processing_flag(table, name):
    res = table.search(Query()[name].exists())
    if res:
        if len(res) > 1:
            raise Exception("expecting only one flag doc for {}".format(name))
        flag_doc = res[0]
    else:
        doc_id = table.insert({name: False})
        flag_doc = table.get(doc_id=doc_id)
    return flag_doc


@contextmanager
def _processing_step(table, name, reprocess=False):
    flag_doc = _get_processing_flag(table, name)
    if not flag_doc[name] or reprocess:
        flag_doc[name] = False
        table.write_back([flag_doc])
        yield True
        flag_doc[name] = True
        table.write_back([flag_doc])
    else:
        yield False


def _get_extension(path):
    return path.split(".")[-1].lower()


def _scan_file(path):
    doc = {}
    doc["path"] = path
    stat = os.stat(path)
    doc["size"] = stat.st_size
    doc["atime"] = datetime.fromtimestamp(stat.st_atime)
    doc["mtime"] = datetime.fromtimestamp(stat.st_mtime)
    doc["ctime"] = datetime.fromtimestamp(stat.st_ctime)
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
    db = get_tinydb(out_path)
    status_table = db.table("status")
    files_table = db.table("files")
    try:
        ### SCAN
        valid_extensions = set(METADATA_READERS.keys())
        if skip_tiff:
            valid_extensions -= set(["tif", "tiff"])
        if skip_nd2:
            valid_extensions.remove("nd2")
        with _processing_step(status_table, "scan", reprocess=rescan) as process:
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
                        if _get_extension(file) in valid_extensions:
                            doc = _scan_file(path)
                            files_table.upsert(doc, Query().path == path)
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
                            if _get_extension(file) in valid_extensions:
                                path = os.path.join(root, file)
                                if not os.path.exists(path):
                                    skipped.append(path)
                                    continue
                                doc = _scan_file(path)
                                files_table.upsert(doc, Query().path == path)
                                num_files_to_process += 1
                print("skipped: {}".format(skipped))
        ### METADATA
        with _processing_step(status_table, "metadata") as process:
            if process:
                skipped = []
                files_to_process = sorted(
                    files_table.search(~Query().metadata.exists()),
                    key=lambda x: x["size"],
                    reverse=True,
                )
                pbar = tqdm_auto(files_to_process)
                for doc in pbar:
                    path = doc["path"]
                    file = os.path.split(path)[-1]
                    extension = _get_extension(file)
                    pbar.set_postfix(
                        file=_abbreviate_filename(file, 20), skipped=len(skipped)
                    )
                    try:
                        doc["metadata"] = METADATA_READERS[extension](path)
                    except KeyboardInterrupt:
                        raise
                    except:
                        skipped.append(path)
                    else:
                        files_table.write_back([doc])
                        gc.collect()
                print("skipped: {}".format(skipped))
    finally:
        db.close()  # flush
