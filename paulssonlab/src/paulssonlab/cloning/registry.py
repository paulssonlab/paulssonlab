import pandas as pd
import os
import re
from cytoolz import dissoc, excepts
from numbers import Integral
from Bio.Seq import Seq
import Bio.Restriction
import pygsheets
from natsort import natsorted, ns
from paulssonlab.api.google import (
    insert_sheet_rows,
    update_sheet_rows,
    get_drive_by_path,
    ensure_drive_folder,
    upload_drive,
    make_drive_folder,
    copy_drive_file,
    copy_drive_folder,
    list_drive,
    clear_sheet,
    FOLDER_MIMETYPE,
    SHEETS_MIMETYPE,
)
from paulssonlab.cloning.workflow import (
    parse_id,
    rename_ids,
    part_entry_to_seq,
    re_digest_part,
    is_bases,
    ID_REGEX,
)
from paulssonlab.api import regex_key
from paulssonlab.api.benchling import upload_sequence
from paulssonlab.cloning.sequence import DsSeqRecord, anneal, pcr, get_seq
from paulssonlab.cloning.commands.semantics import eval_exprs_by_priority
from paulssonlab.cloning.io import (
    value_to_bytes,
    bytes_to_value,
    value_to_mimetype,
    filename_to_mimetype,
    value_to_extension,
)
from paulssonlab.api.util import PROGRESS_BAR

AMBIGUOUS_MIMETYPES = set(["application/octet-stream"])
DEFAULT_LOOKUP_TYPES = ["oligos", "plasmids", "strains", "parts"]
FOLDER_TYPES = ["maps", "sequencing"]
TYPES_WITH_SEQUENCING = ["plasmids", "strains"]
TYPES_WITH_MAPS = ["plasmids"]
TYPES_WITHOUT_IDS = ["parts"]
BENCHLING_SYNC_TYPES = ["parts", "plasmids", "oligos"]
COLLECTION_REGEX = r"^([^_]+)_(\w+)$"
ENABLE_AUTOMATION_FILENAME = "ENABLE_AUTOMATION.txt"


def _name_mapper_for_prefix(old_prefix, new_prefix):
    mapper = lambda name, is_folder: re.sub(
        rf"(\w*){re.escape(old_prefix)}(.*)", rf"\1{new_prefix}\2", name
    )
    return mapper


class GDriveClient(object):
    def __init__(self, registry_key):
        self.registry_key = registry_key

    def __len__(self):
        return len(self.keys())

    def keys(self):
        raise NotImplementedError

    def __delitem__(self, key):
        self[key] = None

    def __iter__(self):
        return iter(self.keys())

    def items(self):
        for key in self:
            yield key, self[key]

    def __getitem__(self, key):
        raise NotImplementedError

    def __setitem__(self, key, value):
        self.local[key] = value

    def update(self, other):
        for key, value in other.items():
            self[key] = value


class SheetClient(GDriveClient):
    def __init__(self, registry_key, registry, id_column=0):
        super().__init__(registry_key)
        self.gdrive_id = registry.gdrive_ids[registry_key[:2]]
        self.id_column = id_column
        if len(registry_key) == 3:
            self.worksheet_title = registry_key[2]
        elif len(registry_key) not in (2, 3):
            raise ValueError(
                "expected key of form (prefix, type) or (prefix, type, worksheet_title)"
            )
        else:
            self.worksheet_title = None
        self.client = registry.sheets_client
        self._spreadsheet = None
        self._worksheet = None
        self._columns = None
        self.rollback()

    def rollback(self):
        self.local = {}
        self._remote = None
        self._remote_index = None

    @property
    def id_column_name(self):
        return self.columns[self.id_column]

    @property
    def spreadsheet(self):
        if self._spreadsheet is None:
            self._spreadsheet = self.client.open_by_key(self.gdrive_id)
        return self._spreadsheet

    @property
    def worksheet(self):
        if self._worksheet is None:
            if self.worksheet_title is not None:
                self._worksheet = self.spreadsheet.worksheet(
                    "title", self.worksheet_title
                )
            else:
                self._worksheet = self.spreadsheet.worksheet()
        return self._worksheet

    @property
    def remote(self):
        if self._remote is None:
            self._download()
        return self._remote

    @property
    def remote_index(self):
        if self._remote_index is None:
            self._download()
        return self._remote_index

    def _download(self):
        df = self.worksheet.get_as_df(value_render=pygsheets.ValueRenderOption.FORMULA)
        # ensure proper ID formatting
        df.iloc[:, self.id_column] = df.iloc[:, self.id_column].str.replace(
            r"\s+", "", regex=True
        )
        self._remote = df
        self._remote_index = pd.Index(df.iloc[:, self.id_column])

    @property
    def columns(self):
        return self.remote.columns

    def _last_id(self):
        ids = sorted(
            [
                parse_id(id_)
                for id_ in self.keys()
                if id_.startswith(self.registry_key[0])
            ]
        )
        if ids:
            last_id = ids[-1]
        else:
            last_id = (self.registry_key[0], 0)
        return last_id

    def next_id(self):
        last_id = self._last_id()
        return f"{last_id[0]}{last_id[1] + 1}"

    def keys(self):
        # don't include rows with empty IDs
        keys = self.local.keys() | set(self.remote_index) - set([""]) - set(
            k for k in self.local.keys() if self.local[k] is None
        )
        return natsorted(keys, alg=ns.IGNORECASE)

    def __contains__(self, key):
        if key.strip() == "":
            return False
        if key in self.local:
            return self.local[key] is not None
        else:
            return key in self.remote_index

    def __getitem__(self, key):
        # don't allow accesing rows with empty IDs
        if key.strip() == "":
            raise KeyError(key)
        if key in self.local:
            row = self.local[key]
            if row is None:
                raise KeyError(key)
        else:
            row = self.remote.iloc[self.remote_index.get_loc(key)].to_dict()
        return {**row, self.id_column_name: key}

    def _validate(self, row):
        if self.id_column_name in row:
            raise ValueError(f"cannot specify ID column: '{self.id_column_name}'")
        unknown_columns = row.keys() - set(self.columns)
        if unknown_columns:
            raise ValueError(f"unknown columns: {unknown_columns}")

    def __setitem__(self, key, row):
        self._validate(row)
        self.local[key] = row

    def append(self, row):
        id_ = self.next_id()
        self[id_] = row
        return id_

    def upsert(self, row, key_columns=[], apply={}):
        """None values in apply is equivalent to the identity function.

        Missing keys in rows are interpreted as empty strings.
        """
        key_columns = list(set(key_columns) | apply.keys())
        if not key_columns:
            raise ValueError(
                "at least one column must be specified in key_columns or apply"
            )
        identity = lambda x: x
        generate_key = lambda row: tuple(
            (apply.get(col, identity) or identity)(row.get(col, ""))
            for col in key_columns
        )
        key_values = generate_key(row)
        matches = [
            id_ for id_, row_ in self.items() if generate_key(row_) == key_values
        ]
        if len(matches) >= 2:
            raise ValueError(
                f"upsert key must be unique, instead found matches: {matches}"
            )
        elif len(matches) == 1:
            id_ = matches[0]
            new_row = dissoc({**self[id_], **row}, self.id_column_name)
            self[id_] = new_row
            return id_
        id_ = self.append(row)
        return id_

    def commit(self, clobber_existing=True, append_only=True, formulae=True):
        """
        Parameters
        ----------
        clobber_existing : bool, optional
            If a local row shares a key with an existing remote row, remote row values will be
            preserved when local row is missing column keys unless clobber_existing is True.
        append_only : bool, optional
            Ensures that when set to False.
        formulae : bool, optional
            If True, any formula values from the first non-header row will be copied (with their
            row number changed) to any new rows.

        """
        rows_to_update = {}
        rows_to_append = []
        # we need to do this to handle IDs without sequential numeric indices (e.g., part names)
        sorter = lambda x: excepts(ValueError, parse_id)(x) or (None, x)
        for key in sorted(self.local.keys(), key=sorter):
            row = self.local[key]
            row_with_id = {**row, self.id_column_name: key}
            if key in self.remote_index:
                rows_to_update[key] = row_with_id
            else:
                rows_to_append.append(row_with_id)
        self._update_rows(rows_to_update)
        self._append_rows(rows_to_append)
        # TODO: we could update self.remote here, but for now let's do the safer thing,
        # which is to fetch everything again from the server
        self.rollback()

    def _update_rows(self, rows):
        rows = [(self.remote_index.get_loc(key), row) for key, row in rows.items()]
        updates = [
            (
                pygsheets.DataRange(
                    (idx, 1), (idx, len(self.columns)), worksheet=self.worksheet
                ),
                [[row.get(col, "") for col in self.columns]],
            )
            for idx, row in rows
        ]
        update_sheet_rows(self.client.sheet.service, updates)

    def _append_rows(self, rows):
        # increment index once for column header, once for 1-indexing
        insert_sheet_rows(
            self.worksheet, len(self.remote) + 2, rows, columns=self.columns
        )


class ItemProxy(object):
    def __init__(self, obj, name):
        self._obj = obj
        self._name = name

    def __contains__(self, key):
        return getattr(self._obj, f"_{self._name}_contains")(key)

    def __delitem__(self, key):
        return getattr(self._obj, f"_{self._name}_delitem")(key)

    def __iter__(self):
        return getattr(self._obj, f"_{self._name}_iter")()

    def items(self):
        return getattr(self._obj, f"_{self._name}_items")()

    def __getitem__(self, key):
        return getattr(self._obj, f"_{self._name}_getitem")(key)

    def __setitem__(self, key, value):
        return getattr(self._obj, f"_{self._name}_setitem")(key, value)

    def update(self, other):
        setitem = getattr(self._obj, f"_{self._name}_setitem")
        for key, value in other.items():
            setitem(key, value)


class FileClient(GDriveClient):
    def __init__(self, registry_key, registry, trash_name="_Trash"):
        super().__init__(registry_key)
        self.gdrive_id = registry.gdrive_ids[registry_key[:2]]
        self.client = registry.drive_service
        self.raw = ItemProxy(self, "raw")
        self.bytes = ItemProxy(self, "bytes")
        self.content = ItemProxy(self, "content")
        self.rollback()
        self._trash_name = trash_name

    def rollback(self):
        self.local = {}
        self._remote = None
        self._remote_folders = None

    @property
    def remote(self):
        if self._remote is None:
            self._download()
        return self._remote

    @property
    def remote_folders(self):
        if self._remote_folders is None:
            self._download()
        return self._remote_folders

    def _download(self):
        self._remote = {}
        self._remote_folders = {}
        self._list_folder((), self.gdrive_id)
        # adding root id as an empty key makes some recursion in commit() easier
        self._remote_folders[()] = {"id": self.gdrive_id}

    def _list_folder(self, key, root):
        files = list_drive(self.client, root=root, fields=["parents"])
        for name, file in files.items():
            new_key = key + (name,)
            if new_key[0] == self._trash_name:
                # don't include trash folder or its contents
                continue
            if file["mimeType"] == FOLDER_MIMETYPE:
                self._remote_folders[new_key] = file
                self._list_folder(new_key, file["id"])
            else:
                self._remote[new_key] = file

    def _download_bytes(self, file):
        if "content" not in file and "bytes" not in file:
            bytes_ = self.client.files().get_media(fileId=file["id"]).execute()
            file["bytes"] = bytes_
            return bytes_
        else:
            return None

    def _get_bytes(self, file):
        if "bytes" not in file:
            if "content" not in file:
                raise ValueError(
                    f"found neither 'bytes' nor 'content' in file: '{file}'"
                )
            else:
                file["bytes"] = value_to_bytes(file["content"])
        return file["bytes"]

    def _get_content(self, file):
        if "content" not in file:
            if "bytes" not in file or "mimeType" not in file:
                raise ValueError(
                    f"found neither 'bytes'/'mimeType' nor 'content' in file: '{file}'"
                )
            else:
                # if mimetype does not indicate file type
                # (e.g., file was uploaded manually using google drive web interface)
                # use filename to infer true mimetype
                if file["mimeType"] in AMBIGUOUS_MIMETYPES:
                    mimetype = filename_to_mimetype(file["name"])
                else:
                    mimetype = file["mimeType"]
                file["content"] = bytes_to_value(file["bytes"], mimetype)
        return file["content"]

    def find(self, key):
        return self._find_loc(key)[0]

    def _find_loc(self, key):
        key = self._validate_key(key)
        for store, is_remote in ((self.local, False), (self.remote, True)):
            filtered_keys = [k[-1] for k in store.keys() if k[:-1] == key[:-1]]
            try:
                matching_key = regex_key(
                    filtered_keys, f"{key[-1]}($|\\.[^\\.]+$)", check_duplicates=True
                )
                full_key = (*key[:-1], matching_key)
                return full_key, is_remote
            except ValueError:
                continue
        raise KeyError(key)

    def _get_loc(self, key):
        full_key, is_remote = self._find_loc(key)
        if not is_remote:
            return self.local[full_key], is_remote
        else:
            return self.remote[full_key], is_remote

    def __contains__(self, key):
        try:
            self._find_loc(key)
            return True
        except KeyError:
            return False

    def _keys_loc(self):
        keys_to_delete = set(k for k, v in self.local.items() if v is None)
        local_keys = set((k, False) for k in self.local.keys())
        remote_keys = set((k, True) for k in self.remote.keys())
        keys = local_keys | remote_keys
        keys = [(k, is_remote) for k, is_remote in keys if k not in keys_to_delete]
        keys = set(sorted(keys, key=lambda k: k[0]))
        return keys

    def keys(self):
        keys = set(k[0] for k in self._keys_loc())
        return natsorted(keys, alg=ns.IGNORECASE)

    def _items_loc(self):
        for key, is_remote in self._keys_loc():
            if is_remote:
                store = self.remote
            else:
                store = self.local
            yield key, store[key], is_remote

    def items(self):
        yield from self._content_items()

    def __getitem__(self, key):
        file, is_remote = self._get_loc(key)
        if is_remote:
            self._download_bytes(file)
        return self._get_content(file)

    def _validate_key(self, key):
        if not isinstance(key, tuple):
            key = (key,)
        if key == ():
            raise KeyError(key)
        return key

    def __setitem__(self, key, value):
        key = self._validate_key(key)
        if isinstance(value, dict):
            for k, v in value.items():
                self[key + (k,)] = v
            return
        try:
            existing_key, is_remote = self._find_loc(key)
        except KeyError:
            existing_key = None
        extension = value_to_extension(value)
        if extension:
            full_key = (*key[:-1], f"{key[-1]}.{extension}")
        else:
            full_key = key
        if value is None:
            self.content[existing_key] = None
        else:
            if existing_key is not None and existing_key != full_key:
                # need to delete existing_key first
                self.content[existing_key] = None
            self.content[full_key] = value

    def _raw_contains(self, key):
        if not isinstance(key, tuple):
            key = (key,)
        if key in self.local:
            return self.local[key] is not None
        else:
            return key in self.remote

    def _raw_delitem(self, key):
        self._raw_setitem(key, None)

    def _raw_iter(self):
        yield from self

    def _raw_items(self):
        for key, value, _ in self._items_loc():
            yield key, value

    def _raw_getitem(self, key):
        return self._raw_getitem_loc(key)[0]

    def _raw_getitem_loc(self, key):
        key = self._validate_key(key)
        if key in self.local:
            value = self.local[key]
            if value is None:
                raise KeyError(key)
            return value, False
        else:
            return self.remote[key], True

    def _raw_setitem(self, key, value):
        key = self._validate_key(key)
        if value is not None:
            if not isinstance(value, dict):
                raise ValueError(f"got {value.__class__}, expecting dict")
            if "mimeType" not in value:
                raise ValueError(f"mimeType is required")
            if not ("content" in value or "bytes" in value):
                raise ValueError(f"at least one of content/bytes is required")
            extra_keys = value.keys() - set(["mimeType", "bytes", "content"])
            if extra_keys:
                raise ValueError(f"got unexpected keys: {extra_keys}")
        base_key = (*key[:-1], os.path.splitext(key[-1])[0])
        try:
            full_key, is_remote = self._find_loc(base_key)
            if full_key != key:
                raise ValueError(
                    f"attempting to set key {key} that conflicts with existing key {full_key} (remote={is_remote})"
                )
        except KeyError:
            pass
        self.local[key] = value

    _bytes_contains = _raw_contains
    _bytes_delitem = _raw_delitem
    _bytes_iter = _raw_iter

    def _bytes_items(self):
        for key, value, is_remote in self._items_loc():
            if is_remote:
                self._download_bytes(value)
            yield key, self._get_bytes(value)

    def _bytes_getitem(self, key):
        file, is_remote = self._raw_getitem_loc(key)
        if is_remote:
            self._download_bytes(file)
        return self._get_bytes(file)

    def _bytes_setitem(self, key, value):
        key = self._validate(key)
        if isinstance(value, dict):
            for k, v in value.items():
                self._bytes_setitem(key + (k,), v)
            return
        if value is None:
            self.local[key] = None
        else:
            if not isinstance(value, bytes):
                raise ValueError("expecting value of type bytes")
            file = self.local.setdefault(key, {})
            file.pop("content", None)
            file["bytes"] = value
            file["mimeType"] = filename_to_mimetype(key[-1])

    _content_contains = _raw_contains
    _content_delitem = _raw_delitem
    _content_iter = _raw_iter

    def _content_items(self):
        for key, value, is_remote in self._items_loc():
            if is_remote:
                self._download_bytes(value)
            yield key, self._get_content(value)

    def _content_getitem(self, key):
        file, is_remote = self._raw_getitem_loc(key)
        if is_remote:
            self._download_bytes(file)
        return self._get_content(file)

    def _content_setitem(self, key, value):
        key = self._validate_key(key)
        if isinstance(value, dict):
            for k, v in value.items():
                self._content_setitem(key + (k,), v)
            return
        if value is None:
            self.local[key] = None
        else:
            file = self.local.setdefault(key, {})
            file.pop("bytes", None)
            file["content"] = value
            file["mimeType"] = value_to_mimetype(value)

    def set_from_file(self, key, path):
        key = self._validate_key(key)
        mimetype = filename_to_mimetype(path)
        with open(path, "rb") as f:
            bytes_ = f.read()
        extension = os.path.splitext(path)[1][1:].lower()
        key = (*key[:-1], f"{key[-1]}.{extension}")
        self.raw[key] = {"bytes": bytes_, "mimeType": mimetype}
        return key

    def update(self, other):
        for key, value in other.items():
            self[key] = value

    def _trash_folder_id(self):
        trash_folder = self.remote_folders.get((self._trash_name,))
        if trash_folder is None:
            trash_folder = make_drive_folder(
                self.client, self._trash_name, self.gdrive_id
            )
            self.remote_folders[(self._trash_name,)] = trash_folder
        return trash_folder["id"]

    def commit(self, overwrite=True, trash="_Trash", remove_empty_folders=True):
        keys_to_trash = set()
        # REMOVE files/folders = None in local
        for key, value in self.local.items():
            if value is None:
                if key in self.remote or key in self.remote_folders:
                    keys_to_trash.add(key)
        for key in keys_to_trash:
            file = self.remote.get(key)
            is_folder = False
            if file is None:
                file = self.remote_folders.get(key)
                is_folder = True
            if file is None:
                continue
            self.client.files().update(
                fileId=file["id"],
                addParents=self._trash_folder_id(),
                removeParents=",".join(file["parents"]),
            ).execute()
            if is_folder:
                # remove all remote files/folders with a key prefixed by "key"
                for k in self.remote:
                    if key == k[: len(key)]:
                        del self.remote[k]
                for k in self.remote_folders:
                    if key == k[: len(key)]:
                        del self.remote_folders[k]
            else:
                del self.remote[key]
        # MAKE FOLDERS for all parent keys for all local keys
        # list() makes a copy so we can delete keys as we iterate
        for key, value in list(self.local.items()):
            if value is not None:
                for i in range(1, len(key)):
                    folder_key = key[:i]
                    if folder_key not in self.remote_folders:
                        folder = make_drive_folder(
                            self.client,
                            folder_key[-1],
                            self.remote_folders[folder_key[:-1]]["id"],
                        )
                        self.remote_folders[folder_key] = folder
                if key in self.remote:
                    if not overwrite:
                        raise ValueError(f"file already exists remotely: {key}")
                    file_id = self.remote[key]["id"]
                else:
                    file_id = None
                upload_drive(
                    self.client,
                    self.bytes[key],
                    key[-1],
                    mimetype=value["mimeType"],
                    file_id=file_id,
                    parent=self.remote_folders[key[:-1]]["id"],
                )
                self.remote[key] = self.local[key]
                del self.local[key]
        # TRASH ALL remaining folders
        if remove_empty_folders:
            nonempty_folders = set([()])  # never remove root
            for key in self.remote:
                for i in range(1, len(key)):
                    nonempty_folders.add(key[:i])
            empty_folders = set(self.remote_folders) - nonempty_folders
            for key in empty_folders:
                folder = self.remote_folders[key]
                self.client.files().update(
                    fileId=folder["id"],
                    addParents=self._trash_folder_id(),
                    removeParents=",".join(folder["parents"]),
                ).execute()
                del self.remote_folders[key]


TYPE_TO_CLIENT = {"maps": FileClient, "sequencing": FileClient, "_default": SheetClient}


class Registry(object):
    def __init__(self, sheets_client, registry_folder, benchling_folder=None):
        self.sheets_client = sheets_client
        self.registry_folder = registry_folder
        self.benchling_folder = benchling_folder
        self.rollback()
        self.refresh()

    @property
    def benchling_session(self):
        return self.benchling_folder.session

    @property
    def sheet_service(self):
        return self.sheets_client.sheet.service

    @property
    def drive_service(self):
        return self.sheets_client.drive.service

    def refresh(self):
        collection_folders = list_drive(
            self.drive_service, root=self.registry_folder, is_folder=True
        )
        gdrive_ids = {}
        for collection_folder in collection_folders.values():
            new_gdrive_ids = self.get_collection(collection_folder["id"])
            duplicate_keys = gdrive_ids.keys() & new_gdrive_ids.keys()
            if len(duplicate_keys):
                raise ValueError(f"found duplicate prefixes: {list(duplicate_keys)}")
            gdrive_ids = {**gdrive_ids, **new_gdrive_ids}
        self.gdrive_ids = gdrive_ids

    def get_collection(self, collection_folder):
        files = list_drive(self.drive_service, root=collection_folder)
        if ENABLE_AUTOMATION_FILENAME not in files:
            return {}
        gdrive_ids = {}
        for file in files.values():
            name = file["name"]
            match = re.match(COLLECTION_REGEX, name)
            if match:
                prefix = match.group(1)
                type_ = match.group(2)
                ensure_drive_folder(file, type_ in FOLDER_TYPES)
                key = (prefix, type_)
                if key in gdrive_ids:
                    raise ValueError(f"found duplicate prefix: {key}")
                gdrive_ids[key] = file["id"]
        return gdrive_ids

    def __iter__(self):
        yield from self.keys()

    def keys(self):
        return self.gdrive_ids.keys()

    def items(self):
        for key in self:
            yield key, self[key]

    def __contains__(self, key):
        return key in self.keys()

    def __getitem__(self, key):
        if key in self.clients:
            return self.clients[key]
        else:
            client_type = TYPE_TO_CLIENT.get(key[1], TYPE_TO_CLIENT["_default"])
            client = client_type(key, self)
            self.clients[key] = client
            return client

    def update(self, other):
        for key, value in other.items():
            self[key].update(value)

    def commit(self):
        for client in self.clients.values():
            client.commit()

    def rollback(self):
        self.clients = {}

    def duplicate_collection(
        self,
        source_prefix,
        dest_prefix,
        source_folder_name=None,
        dest_folder_name=None,
        clear=True,
    ):
        # TODO: handle parts spreadsheet clearing (keep formulae)
        if source_folder_name is None:
            source_folder_name = f"{source_prefix}_collection"
        if dest_folder_name is None:
            dest_folder_name = f"{dest_prefix}_collection"
        collections = list_drive(self.drive_service, root=self.registry_folder)
        if dest_folder_name in collections:
            raise ValueError(f"collection '{dest_folder_name}' already exists")
        if source_folder_name not in collections:
            raise ValueError(f"collection '{source_folder_name}' not found")
        source_folder = collections[source_folder_name].get("id")
        source_files = list_drive(self.drive_service, root=source_folder)
        dest_folder = make_drive_folder(
            self.drive_service, dest_folder_name, self.registry_folder
        )["id"]
        for source_file in source_files.values():
            if source_file["mimeType"] == FOLDER_MIMETYPE:
                continue
            dest_file_name = None
            if source_file["name"] == ENABLE_AUTOMATION_FILENAME:
                dest_file_name = ENABLE_AUTOMATION_FILENAME
            else:
                match = re.match(COLLECTION_REGEX, source_file["name"])
                if match:
                    dest_type_prefix = re.sub(
                        f"{re.escape(source_prefix)}$", dest_prefix, match.group(1)
                    )
                    dest_file_name = f"{dest_type_prefix}_{match.group(2)}"
            if dest_file_name is not None:
                dest_body = {"name": dest_file_name, "parents": [dest_folder]}
                dest_file = (
                    self.drive_service.files()
                    .copy(fileId=source_file["id"], body=dest_body)
                    .execute()
                )
                if match:
                    source_type_prefix = match.group(1)
                    source_type = match.group(2)
                    name_mapper = _name_mapper_for_prefix(
                        source_type_prefix, dest_type_prefix
                    )
                    if source_type in TYPES_WITH_SEQUENCING:
                        source_seq_folder_name = f"{source_type_prefix}_sequencing"
                        dest_seq_folder_name = f"{dest_type_prefix}_sequencing"
                        dest_seq_folder = make_drive_folder(
                            self.drive_service, dest_seq_folder_name, dest_folder
                        )["id"]
                        if source_seq_folder_name in source_files and not clear:
                            copy_drive_folder(
                                self.drive_service,
                                source_files[source_seq_folder_name]["id"],
                                dest_seq_folder,
                                transform_names=name_mapper,
                            )
                    if source_type in TYPES_WITH_MAPS:
                        source_map_folder_name = f"{source_type_prefix}_maps"
                        dest_map_folder_name = f"{dest_type_prefix}_maps"
                        dest_map_folder = make_drive_folder(
                            self.drive_service, dest_map_folder_name, dest_folder
                        )["id"]
                        if source_map_folder_name in source_files and not clear:
                            copy_drive_folder(
                                self.drive_service,
                                source_files[source_map_folder_name]["id"],
                                dest_map_folder,
                                transform_names=name_mapper,
                            )
                    if source_file["mimeType"] == SHEETS_MIMETYPE:
                        if clear:
                            # get first worksheet
                            dest_sheet = self.sheets_client.open_by_key(
                                dest_file["id"]
                            ).worksheet()
                            clear_sheet(dest_sheet)
                        elif source_type not in TYPES_WITHOUT_IDS:
                            # get first worksheet
                            dest_sheet = self.sheets_client.open_by_key(
                                dest_file["id"]
                            ).worksheet()
                            rename_ids(dest_sheet, source_prefix, dest_type_prefix)
        return dest_folder

    def get_prefix_and_types(self, name):
        match = re.match(ID_REGEX, name)
        if match:
            prefix = match.group(1)
            types = [k[1] for k in self.keys() if k[0] == prefix]
        else:
            prefix = None
            types = ["parts"]
        return prefix, types

    def get_with_prefix(self, name, prefix, type_):
        if type_ == "parts":
            client = None
            for key in self.keys():
                if key[1] != "parts":
                    continue
                if name in self[key]:
                    prefix = key[0]
                    client = self[key]
                    break
            if client is None:
                raise ValueError(f"could not find part '{name}'")
        else:
            key = (prefix, type_)
            if key not in self:
                raise ValueError(f"prefix '{prefix}' does not have a '{type_}'")
            client = self[key]
            if name not in client:
                raise ValueError(f"could not find '{name}'")
        entry = client[name]
        entry["_id"] = name
        entry["_type"] = type_
        entry["_prefix"] = prefix
        if type_ in TYPES_WITH_MAPS:
            try:
                map_client = self[(prefix, "maps")]
            except:
                raise ValueError(f"no maps folder for prefix '{prefix}'")
            if name not in map_client:
                raise ValueError(f"could not find map for '{name}'")
            entry["_seq"] = map_client[name]
        elif type_ == "parts":
            res = eval_exprs_by_priority(entry["Usage"], self.get)
            if res is None:
                seq = None
            else:
                seq = res.get("_seq")
            if seq is None:
                seq = part_entry_to_seq(entry)
            seq.id = seq.name = seq.description = name
            seq = seq.annotate(name)
            entry["_seq"] = seq
        elif "Sequence" in entry:
            # oligos are sometimes ds-DNA and sometimes ss-DNA
            # so make everything into a DsSeqRecord
            # this also enables nice convenience functions like .annotate()
            entry["_seq"] = DsSeqRecord(Seq(entry["Sequence"]))
        return entry

    def get(self, name, types=("plasmids", "strains", "oligos", "parts")):
        match = re.match(ID_REGEX, name)
        if match:
            prefix = match.group(1)
        else:
            prefix = None
            types = ("parts",)
        for type_ in types:
            try:
                return self.get_with_prefix(name, prefix, type_)
            except ValueError:
                pass
        raise ValueError(f"could not find '{name}'")

    def sync_benchling(
        self, overwrite=None, return_data=False, progress_bar=PROGRESS_BAR
    ):
        raise NotImplementedError  # TODO: update this to work with new registry
        cache = {}
        problems = []
        prefixes = self.registry.keys()
        if progress_bar is not None and len(prefixes):
            prefixes = progress_bar(prefixes)
        for prefix, type_ in prefixes:
            if type_ not in BENCHLING_SYNC_TYPES:
                continue
            if hasattr(prefixes, "set_description"):
                prefixes.set_description(f"Scanning {prefix} ({type_})")
            oligo = type_ == "oligos"
            df = self.get_df((prefix, type_))
            rows = df.index[:1]  # TODO
            if progress_bar is not None and len(rows) >= 2:
                rows = progress_bar(rows)
            for name in rows:
                if not name:
                    continue
                if hasattr(rows, "set_description"):
                    rows.set_description(f"Processing {name}")
                try:
                    seq = self._get(df, prefix, type_, name)["_seq"]
                except Exception as e:
                    problems.append((prefix, type_, name, e))
                    continue
                bases = str(get_seq(seq))
                if not is_bases(bases):
                    problems.append((prefix, type_, name, None))
                    continue
                # for parts, add annotations to indicate overhangs
                # because Benchling can't display dsDNA with sticky ends
                if type_ == "parts":
                    seq = seq.annotate_overhangs()
                path = (f"{prefix}_{type_}", name)
                dna = upload_sequence(
                    self.benchling_folder,
                    path,
                    seq,
                    oligo=oligo,
                    overwrite=overwrite,
                    cache=cache,
                )
                if return_data:
                    cache[path] = dna
        if return_data:
            return problems, cache
        else:
            return problems
