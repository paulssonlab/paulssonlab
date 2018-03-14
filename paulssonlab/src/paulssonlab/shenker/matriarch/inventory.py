import numpy as np
import click
from tinydb import TinyDB
from datetime import datetime
from collections import OrderedDict
import os
from util import tqdm_auto
from metadata import read_nd2_file_metadata

read_tiff_file_metadata = lambda p: None
METADATA_READERS = {
    "tif": read_tiff_file_metadata,
    "tiff": read_tiff_file_metadata,
    "nd2": read_nd2_file_metadata,
}


@click.command()
@click.argument("in_path", type=click.Path(exists=True))
@click.argument("out_path", type=click.Path(writable=True, dir_okay=False))
@click.option("--metadata/--no-metadata", default=True)
def inventory(in_path, out_path, metadata):
    db = TinyDB(out_path)
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
