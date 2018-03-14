import numpy as np
import click
from tinydb import TinyDB
from tinydb_serialization import SerializationMiddleware, Serializer
import pickle
import base64
import zarr
from datetime import datetime
from collections import OrderedDict
import os
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


@click.command()
@click.argument("in_path", type=click.Path(exists=True))
@click.argument("out_path", type=click.Path(writable=True, dir_okay=False))
@click.option("--metadata/--no-metadata", default=True)
def inventory(in_path, out_path, metadata):
    serialization = SerializationMiddleware()
    # serialization.register_serializer(DateTimeSerializer(), 'datetime')
    serialization.register_serializer(ZarrSerializer(), "zarr")
    serialization.register_serializer(BytesSerializer(), "bytes")
    db = TinyDB(out_path, storage=serialization)
    files_to_process = OrderedDict()
    valid_extensions = set(METADATA_READERS.keys())
    for root, dirs, files in tqdm_auto(
        os.walk(in_path), desc="scanning for image files"
    ):
        # modifying dirs in place will affect which dirs are traversed
        # SEE: https://stackoverflow.com/questions/19859840/excluding-directories-in-os-walk
        dirs[:] = [
            d for d in dirs if not os.path.exists(os.path.join(root, d, ".zattrs"))
        ]
        for file in files:
            extension = file.split(".")[-1].lower()
            if extension in valid_extensions:
                path = os.path.join(root, file)
                files_to_process[path] = METADATA_READERS[extension]
    for path, metadata_reader in tqdm_auto(files_to_process.items()):
        entry = {}
        entry["path"] = path
        stat = os.stat(path)
        entry["size"] = stat.st_size
        entry["atime"] = datetime.fromtimestamp(stat.st_atime)
        entry["mtime"] = datetime.fromtimestamp(stat.st_mtime)
        entry["ctime"] = datetime.fromtimestamp(stat.st_ctime)
        if metadata:
            entry["metadata"] = metadata_reader(path)
        db.insert(entry)
