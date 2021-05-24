import pandas as pd
import re
from numbers import Integral
from Bio.Seq import Seq
import Bio.Restriction
import pygsheets
from paulssonlab.api.google import (
    insert_sheet_rows,
    update_sheet_rows,
    get_drive_by_path,
    ensure_drive_folder,
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
)
from paulssonlab.api import regex_key
from paulssonlab.api.benchling import upload_sequence
from paulssonlab.cloning.sequence import DsSeqRecord, anneal, pcr, get_seq
from paulssonlab.cloning.commands.semantics import eval_exprs_by_priority
from paulssonlab.cloning.io import read_sequence
from paulssonlab.api.util import PROGRESS_BAR

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

    def keys(self):
        raise NotImplementedError

    def items(self):
        for key in self.keys():
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
        self.clear_cache()

    def clear_cache(self):
        self.local = {}
        self._remote = None
        self._remote_index = None

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
        return self.local.keys() | set(self.remote_index) - set([""])

    def __getitem__(self, key):
        # don't allow accesing rows with empty IDs
        if key.strip() == "":
            raise KeyError(key)
        if key in self.local:
            row = self.local[key]
        else:
            row = self.remote.iloc[self.remote_index.get_loc(key)].to_dict()
        return {**row, self.columns[self.id_column]: key}

    def _validate(self, row):
        if self.columns[self.id_column] in row:
            raise ValueError(
                f"cannot specify ID column: '{self.columns[self.id_column]}'"
            )
        unknown_columns = row.keys() - set(self.columns)
        if unknown_columns:
            raise ValueError(f"unknown columns: {unknown_columns}")

    def __setitem__(self, key, row):
        # ensure ID is parsable
        parse_id(key)
        self._validate(row)
        self.local[key] = row

    def append(self, row):
        id_ = self.next_id()
        self[id_] = row
        return id_

    def upsert(self, key_columns, row):
        key_values = tuple(row[col] for col in key_columns)
        key = None
        for location in ("local", "remote"):
            store = getattr(self, location)
            if location == "remote":
                rows = df.iterrows()
            else:
                rows = store.items()
            matches = [
                key
                for key, row in rows
                if tuple(row[col] for col in key_columns) == key_values
            ]
            if len(matches) >= 2:
                raise ValueError(
                    f"upsert key must be unique, instead found {location} matches: {matches}"
                )
            elif len(matches) == 1:
                key = matches[0]
                self[key] = {**store[key], **row}
                return key
        key = self.next_id()
        self.local[key] = row
        return key

    def save(self, clobber_existing=True, append_only=True, formulae=True):
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
        id_column = self.columns[self.id_column]
        for key in sorted(self.local.keys(), key=parse_id):
            row = self.local[key]
            row_with_id = {**row, id_column: key}
            if key in self.remote_index:
                rows_to_update[key] = row_with_id
            else:
                rows_to_append.append(row_with_id)
        self._update_rows(rows_to_update)
        self._append_rows(rows_to_append)
        # TODO: we could update self.remote here, but for now let's do the safer thing,
        # which is to fetch everything again from the server
        self.clear_cache()

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

    def __getitem__(self, key):
        return getattr(self._obj, f"_{self._name}_getitem")(self._obj, key)

    def __setitem__(self, key, value):
        return getattr(self._obj, f"_{self._name}_setitem")(self._obj, key, value)

    def update(self, other):
        setitem = getattr(self._obj, f"_{self._name}_setitem")
        for key, value in other.items():
            self[key] = value


class FileClient(GDriveClient):
    def __init__(self, registry_key, registry):
        super().__init__(registry_key)
        self.gdrive_id = registry.gdrive_ids[registry_key[:2]]
        self.client = registry.drive_service
        self.raw = ItemProxy(self, "raw")
        self.bytes = ItemProxy(self, "bytes")
        self.content = ItemProxy(self, "content")
        self.clear_cache()

    def clear_cache(self):
        self._remote = None

    @property
    def remote(self):
        if self._remote is None:
            self._download_file_list()
        return self._remote

    def _download_file_list(self):
        self._remote = list_drive(self.client, root=self.gdrive_id)

    def _download_bytes(self, file):
        bytes_ = self.client.files().get_media(fileId=file["id"]).execute()
        file["bytes"] = bytes_
        return bytes_
        # if file["kind"] == "drive#folder":
        #     files = list_drive(self.client, root=file["id"])
        #     # folder_content = {name: self._get_file(f) for name, f in files.items()}
        #     return folder_content
        # else:
        #     content = (
        #         self.client.files()
        #         .get_media(fileId=file["id"])
        #         .execute()
        #     )
        # try:
        #     content = read_sequence(content, filename=file["name"])
        # except:
        #     pass
        # return content

    def keys(self):
        return self.local.keys() | self.remote.keys()

    def key(self, key):
        return self._find(key)[0]

    def _find(self, key):
        for store in (self.local, self.remote):
            try:
                return regex_key(
                    store, f"{key}($|\\.)", check_duplicates=True, return_key=True
                )
            except:
                pass
        raise KeyError(key)

    def _get_bytes(self, file):
        if "bytes" not in file:
            if "content" not in file:
                self._download_bytes(file)
            else:
                file["bytes"] = object_to_bytes(file["content"])
        return file["bytes"]

    def _get_content(self, file):
        if "content" not in file:
            self._download_bytes(file)
            file["content"] = bytes_to_object(file["bytes"])
        return file["content"]

    def __getitem__(self, key):
        return self._get_content(self._find(key)[1])

    def __setitem__(self, key, value):
        extension = object_to_extension(value)
        full_key = f"{key}.{extension}"
        self.content[full_key] = value

    def _raw_getitem(self, key):
        if key in self.local:
            return self.local[key]
        else:
            return self.remote[key]

    def _raw_setitem(self, key, value):
        self.local[key] = value
        # TODO: validate value is dict, contains
        # - at least one of content/bytes
        # - mimetype
        # - silently delete id?
        # - no other keys

    def _bytes_getitem(self, key):
        return self._get_bytes(self.raw[key])
        # value = self._find(key)[1]
        # if "bytes" in value:
        #     return value["bytes"]
        # elif "content" in value:
        #     return to_bytes(value["content"])
        # else:
        # must download if remote?
        # to_bytes()

    def _bytes_setitem(self, key, value):
        file = self.local.setdefault(key, {})
        file.pop("content", None)
        file["bytes"] = value
        if not isinstance(value, bytes):
            raise ValueError("expecting object of type bytes")
        file["mimetype"] = filename_to_mimetype(key)

    def _content_getitem(self, key):
        return self._get_content(self.raw[key])

    def _content_setitem(self, key, value):
        file = self.local.setdefault(key, {})
        file.pop("bytes", None)
        file["content"] = value
        file["mimetype"] = object_to_mimetype(value)

    def set_from_file(self, path):
        pass
        # basename

    def update(self, other):
        for key, value in other.items():
            self[key] = value

    def save(self, overwrite=True):
        pass
        # for key, content in self.local.items():
        #     try:
        #         existing_file = self._find_file(key)
        #     except:
        #         existing_file = None
        #     if overwrite and existing_file is not None:
        #         pass
        #         # delete file
        #     if content is not None:
        #         pass

        # upload
        # None value in local means delete
        # delete name\..* before uploading name.gb
        # add .gb (k?) extension for DsSeqRecord/SeqRecord, .fasta for Seq (?)
        # TODO: pass through raw content/mimetype/etc.? wrap as GDriveFile/dict for folder?

        # def clear_cache(self):
        #     self._remote = None
        # self._files = None
        # self.remote = {}

        # @property
        # def files(self):
        #     if self._files is None:
        #         self._files = list_drive(self.client, root=self.gdrive_id)
        #     return self._files

        # def __getitem__(self, key):
        #     if key in self.local:
        #         return self.local[key]
        #     elif key in self.remote:
        #         return self.remote[key]
        #     else:
        #         file = self._find_file(key)
        #         if file is None:
        #             raise ValueError(f"key '{key}' not found")
        #         content = self._fetch_file(file)
        #         self.remote[key] = content
        #         return content

        # def _find_file(self, key):
        #     return regex_key(self.files, f"{key}($|\\.)", check_duplicates=True)

        # def _fetch_file(self, file):
        #     if file["kind"] == "drive#folder":
        #         files = list_drive(self.client, root=file["id"])
        #         folder_content = {name: self._get_file(f) for name, f in files.items()}
        #         return folder_content
        #     else:
        #         content = (
        #             self.client.files()
        #             .get_media(fileId=file["id"])
        #             .execute()
        #             .decode("utf8")
        #         )
        #         try:
        #             content = read_sequence(content, filename=file["name"])
        #         except:
        #             pass
        #         return content

        # def update(self, other):
        #     for key, value in other.items():
        #         self[key] = value

        # def save(self, overwrite=True):
        #     for key, content in self.local.items():
        #         try:
        #             existing_file = self._find_file(key)
        #         except:
        #             existing_file = None
        #         if overwrite and existing_file is not None:
        #             pass
        #             # delete file
        #         if content is not None:
        #             pass
        # upload
        # None value in local means delete
        # delete name\..* before uploading name.gb
        # add .gb (k?) extension for DsSeqRecord/SeqRecord, .fasta for Seq (?)
        # TODO: pass through raw content/mimetype/etc.? wrap as GDriveFile/dict for folder?


TYPE_TO_CLIENT = {"maps": FileClient, "sequencing": FileClient, "_default": SheetClient}


class Registry(object):
    def __init__(self, sheets_client, registry_folder, benchling_folder=None):
        self.sheets_client = sheets_client
        self.registry_folder = registry_folder
        self.benchling_folder = benchling_folder
        self.clear_cache()
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

    @property
    def registry(self):
        return self.gdrive_ids.keys()

    def __iter__(self):
        yield from self.registry

    def items(self):
        for key in self:
            yield key, self[key]

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

    def save(self):
        for client in self.clients.values():
            client.save()

    def clear_cache(self):
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
        )
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
                        )
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
                        )
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

    # def get_sheet_by_id(self, id_):
    #     sheet = self.sheets.get(id_)
    #     if sheet is None:
    #         sheet = self.sheets_client.open_by_key(id_).worksheet()
    #         self.sheets[id_] = sheet
    #     return sheet

    # def get_sheet_by_id(self, id_):
    #         sheet = self.sheets.get(id_)
    #         if sheet is None:
    #             if len(id_) not in (2, 3):
    #                 raise ValueError(
    #                     "expecting id tuple of form (prefix, type) or (prefix, type, sheet_index/sheet_title)"
    #                 )
    #             prefix_type = id_[:2]
    #             if len(id_) == 3:
    #                 sheet_id = id_[2]
    #             else:
    #                 sheet_id = None
    #             if sheet_id is None:
    #                 sheet_id = 0
    #             if isinstance(sheet_id, Integral):
    #                 sheet_id_type = "index"
    #             else:
    #                 sheet_id_type = "title"
    #             sheet = self.sheets_client.open_by_key(prefix_type).worksheet(
    #                 sheet_id_type, sheet_id
    #             )
    #             self.sheets[id_] = sheet
    #         return sheet

    # def get_sheet(self, key):
    #     return self.get_sheet_by_id(self.registry[key])

    # def get_df_by_id(self, id_):
    #     df = self.dfs.get(id_)
    #     if df is None:
    #         df = self.get_sheet_by_id(id_).get_as_df(index_column=1)
    #         self.dfs[id_] = df
    #     return df

    # def get_df_by_name(self, name, types=("plasmids", "strains", "oligos", "parts")):
    #     match = re.match(ID_REGEX, name)
    #     if match:
    #         prefix = match.group(1)
    #     else:
    #         types = ("parts",)
    #     for type_ in types:
    #         if type_ == "parts":
    #             for (prefix, part_type), id_ in self.registry.items():
    #                 if part_type != "parts":
    #                     continue
    #                 sheet = self.get_df_by_id(id_)
    #                 try:
    #                     idx = sheet.index.get_loc(name)
    #                     return sheet, prefix, part_type
    #                 except KeyError:
    #                     continue
    #             break
    #         else:
    #             id_ = self.registry.get((prefix, type_))
    #             if id_ is not None:
    #                 return self.get_df_by_id(id_), prefix, type_
    #     raise ValueError(f"could not find dataframe for '{name}'")

    # def get_df(self, key):
    #     return self.get_df_by_id(self.registry[key])

    # def get_next_empty_row(worksheet, skip_columns=0):
    #     last_idx, _ = _get_next_empty_row(worksheet, skip_columns=skip_columns)
    #     if last_idx is None:
    #         return 2
    #     else:
    #         # increment twice for:
    #         # - add one to convert from zero-indexing to one-indexing
    #         # - row 1 is header
    #         return last_idx + 2

    # def get_empty_column_mask(key):
    #     if key in self.column_masks:
    #         return self.column_masks[key]
    #     else:
    #         worksheet = self.get_sheet(key)
    #         mask = empty_column_mask(worksheet)
    #         self.column_masks[key] = mask
    #         return mask

    # def _get_next_empty_row(worksheet, skip_columns=0):
    #     df = worksheet.get_as_df(value_render=pygsheets.ValueRenderOption.FORMULA)
    #     formula_mask = df.iloc[0].str.startswith("=").values
    #     has_datavalidation = columns_with_validation(
    #         worksheet.client.sheet.service,
    #         worksheet.spreadsheet.id,
    #         worksheet.spreadsheet._sheet_list,
    #     )
    #     validation_mask = has_datavalidation[worksheet.title]
    #     validation_mask += [False] * (len(df.columns) - len(validation_mask))
    #     validation_mask = np.array(validation_mask)
    #     mask = formula_mask | validation_mask
    #     masked_values = df.iloc[:, skip_columns:].iloc[:, ~mask[skip_columns:]]
    #     masked_values[masked_values == ""] = np.nan
    #     nonempty = ~masked_values.isnull().all(axis=1)
    #     last_idx = nonempty[nonempty].last_valid_index()
    #     if last_idx is not None:
    #         # convert to Python int because
    #         # DataFrame.last_valid_index() returns np.int64, which is not JSON-serializable
    #         last_idx = int(last_idx)
    #     return last_idx, df

    # def get_next_collection_id(worksheet):
    #     last_idx, df = _get_next_empty_row(worksheet, skip_columns=1)
    #     prefix = worksheet.spreadsheet.title.split("_")[0]
    #     if last_idx is None:
    #         num = 0
    #     else:
    #         num = last_idx + 1
    #     # increment twice for:
    #     # - add one to convert from zero-indexing to one-indexing
    #     # - row 1 is header
    #     row = num + 2
    #     return (prefix, num + 1), row

    # def get_next_id(self, key):
    #     if key in self.ids:
    #         id_ = self.ids[key]
    #         (prefix, num), row = id_
    #         self.ids[key] = ((prefix, num + 1), row + 1)
    #         return id_
    #     sheet = self.get_sheet(key)

    # def types_for_prefix(self, prefix):
    #     types = []
    #     for reg_prefix, type_ in self.registry.keys():
    #         if prefix == reg_prefix:
    #             types.append(type_)
    #     return types

    # def get_map_list(self, prefix):
    #     if prefix in self.maps:
    #         return self.maps[prefix]
    #     else:
    #         key = (prefix, "maps")
    #         maps_folder = self.registry.get(key)
    #         if maps_folder is None:
    #             raise ValueError(f"expecting a registry entry {key}")
    #         maps = list_drive(
    #             self.drive_service, root=maps_folder, is_folder=False
    #         )
    #         self.maps[prefix] = maps
    #         return maps

    # def _get_map(self, prefix, name):
    #     seq_files = self.get_map_list(prefix)
    #     seq_file = regex_key(seq_files, f"{name}\\.", check_duplicates=True)
    #     seq = read_sequence(
    #         self.drive_service.files()
    #         .get_media(fileId=seq_file["id"])
    #         .execute()
    #         .decode("utf8")
    #     )
    #     molecule_type = seq.annotations["molecule_type"]
    #     if molecule_type == "ds-DNA":
    #         seq = DsSeqRecord(seq, id=name, name=name, description=name)
    #     else:
    #         raise NotImplementedError(
    #             f"cannot import genbank molecule_type: '{molecule_type}'"
    #         )
    #     entry = {"_seq": seq}
    #     return entry

    # def _get(self, df, prefix, type_, name):
    #     if name not in df.index:
    #         raise ValueError(f"cannot find {name}")
    #     entry = df.loc[name].to_dict()
    #     entry["_id"] = name
    #     entry["_type"] = type_
    #     entry["_prefix"] = prefix
    #     if type_ in TYPES_WITH_MAPS:
    #         entry.update(self._get_map(prefix, name))
    #     elif type_ == "parts":
    #         res = eval_exprs_by_priority(entry["Usage"], self.get)
    #         if res is None:
    #             seq = None
    #         else:
    #             seq = res.get("_seq")
    #         if seq is None:
    #             seq = part_entry_to_seq(entry)
    #         seq.id = seq.name = seq.description = name
    #         entry["_seq"] = seq
    #     elif "Sequence" in entry:
    #         entry["_seq"] = Seq(entry["Sequence"])
    #     return entry

    # def get(self, name, types=("plasmids", "strains", "oligos", "parts"), seq=True):
    #     df, prefix, type_ = self.get_df_by_name(name, types=types)
    #     return self._get(df, prefix, type_, name)

    def sync_benchling(
        self, overwrite=None, return_data=False, progress_bar=PROGRESS_BAR
    ):
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
