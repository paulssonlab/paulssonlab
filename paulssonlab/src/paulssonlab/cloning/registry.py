import re
from numbers import Integral
from Bio.Seq import Seq
import Bio.Restriction
from paulssonlab.api.google import (
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
    rename_ids,
    ID_REGEX,
    part_entry_to_seq,
    re_digest_part,
    is_bases,
)
from paulssonlab.api import read_sequence, regex_key
from paulssonlab.api.benchling import upload_sequence
from paulssonlab.cloning.sequence import DsSeqRecord, anneal, pcr, get_seq
from paulssonlab.cloning.commands.semantics import eval_exprs_by_priority
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


def apply_expr(func, *args):
    # this is a placeholder to wrap all sequence expressions as {"_seq": DsSeqRecord(...)}
    # until we have a need to pass information beyond the bare sequence
    args = [a["_seq"] if isinstance(a, dict) and "_seq" in a else a for a in args]
    return {"_seq": func(*args)}


class Registry(object):
    def __init__(self, sheets_client, registry_folder, benchling_folder=None):
        self.clear_cache()
        self.sheets_client = sheets_client
        self.registry_folder = registry_folder
        self.benchling_folder = benchling_folder
        self.refresh()

    @property
    def benchling_session(self):
        return self.benchling_folder.session

    def clear_cache(self):
        self.sheets = {}
        self.dfs = {}
        self.maps = {}
        self.ids = {}

    def refresh(self):
        collection_folders = list_drive(
            self.sheets_client.drive.service, root=self.registry_folder, is_folder=True
        )
        registry = {}
        for collection_folder in collection_folders.values():
            new_registry = self.get_collection(collection_folder["id"])
            duplicate_keys = registry.keys() & new_registry.keys()
            if len(duplicate_keys):
                raise ValueError(f"found duplicate prefixes: {list(duplicate_keys)}")
            registry = {**registry, **new_registry}
        self.registry = registry

    def get_collection(self, collection_folder):
        files = list_drive(self.sheets_client.drive.service, root=collection_folder)
        if ENABLE_AUTOMATION_FILENAME not in files:
            return {}
        registry = {}
        for file in files.values():
            name = file["name"]
            match = re.match(COLLECTION_REGEX, name)
            if match:
                prefix = match.group(1)
                type_ = match.group(2)
                ensure_drive_folder(file, type_ in FOLDER_TYPES)
                key = (prefix, type_)
                if key in registry:
                    raise ValueError(f"found duplicate prefix: {key}")
                registry[key] = file["id"]
        return registry

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
        drive_service = self.sheets_client.drive.service
        collections = list_drive(drive_service, root=self.registry_folder)
        if dest_folder_name in collections:
            raise ValueError(f"collection '{dest_folder_name}' already exists")
        if source_folder_name not in collections:
            raise ValueError(f"collection '{source_folder_name}' not found")
        source_folder = collections[source_folder_name].get("id")
        source_files = list_drive(drive_service, root=source_folder)
        dest_folder = make_drive_folder(
            drive_service, dest_folder_name, self.registry_folder
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
                    drive_service.files()
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
                            drive_service, dest_seq_folder_name, dest_folder
                        )
                        if source_seq_folder_name in source_files and not clear:
                            copy_drive_folder(
                                drive_service,
                                source_files[source_seq_folder_name]["id"],
                                dest_seq_folder,
                                transform_names=name_mapper,
                            )
                    if source_type in TYPES_WITH_MAPS:
                        source_map_folder_name = f"{source_type_prefix}_maps"
                        dest_map_folder_name = f"{dest_type_prefix}_maps"
                        dest_map_folder = make_drive_folder(
                            drive_service, dest_map_folder_name, dest_folder
                        )
                        if source_map_folder_name in source_files and not clear:
                            copy_drive_folder(
                                drive_service,
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

    def get_sheet_by_id(self, id_):
        sheet = self.sheets.get(id_)
        if sheet is None:
            sheet = self.sheets_client.open_by_key(id_).worksheet()
            self.sheets[id_] = sheet
        return sheet

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

    def get_sheet(self, key):
        return self.get_sheet_by_id(self.registry[key])

    def get_df_by_id(self, id_):
        df = self.dfs.get(id_)
        if df is None:
            df = self.get_sheet_by_id(id_).get_as_df(index_column=1)
            self.dfs[id_] = df
        return df

    def get_df_by_name(self, name, types=("plasmids", "strains", "oligos", "parts")):
        match = re.match(ID_REGEX, name)
        if match:
            prefix = match.group(1)
        else:
            types = ("parts",)
        for type_ in types:
            if type_ == "parts":
                for (prefix, part_type), id_ in self.registry.items():
                    if part_type != "parts":
                        continue
                    sheet = self.get_df_by_id(id_)
                    try:
                        idx = sheet.index.get_loc(name)
                        return sheet, prefix, part_type
                    except KeyError:
                        continue
                break
            else:
                id_ = self.registry.get((prefix, type_))
                if id_ is not None:
                    return self.get_df_by_id(id_), prefix, type_
        raise ValueError(f"could not find dataframe for '{name}'")

    def get_df(self, key):
        return self.get_df_by_id(self.registry[key])

    def get_next_empty_row(worksheet, skip_columns=0):
        last_idx, _ = _get_next_empty_row(worksheet, skip_columns=skip_columns)
        if last_idx is None:
            return 2
        else:
            # increment twice for:
            # - add one to convert from zero-indexing to one-indexing
            # - row 1 is header
            return last_idx + 2

    def get_empty_column_mask(key):
        if key in self.column_masks:
            return self.column_masks[key]
        else:
            worksheet = self.get_sheet(key)
            mask = empty_column_mask(worksheet)
            self.column_masks[key] = mask
            return mask

    def _get_next_empty_row(worksheet, skip_columns=0):
        df = worksheet.get_as_df(value_render=pygsheets.ValueRenderOption.FORMULA)
        formula_mask = df.iloc[0].str.startswith("=").values
        has_datavalidation = columns_with_validation(
            worksheet.client.sheet.service,
            worksheet.spreadsheet.id,
            worksheet.spreadsheet._sheet_list,
        )
        validation_mask = has_datavalidation[worksheet.title]
        validation_mask += [False] * (len(df.columns) - len(validation_mask))
        validation_mask = np.array(validation_mask)
        mask = formula_mask | validation_mask
        masked_values = df.iloc[:, skip_columns:].iloc[:, ~mask[skip_columns:]]
        masked_values[masked_values == ""] = np.nan
        nonempty = ~masked_values.isnull().all(axis=1)
        last_idx = nonempty[nonempty].last_valid_index()
        if last_idx is not None:
            # convert to Python int because
            # DataFrame.last_valid_index() returns np.int64, which is not JSON-serializable
            last_idx = int(last_idx)
        return last_idx, df

    def get_next_collection_id(worksheet):
        last_idx, df = _get_next_empty_row(worksheet, skip_columns=1)
        prefix = worksheet.spreadsheet.title.split("_")[0]
        if last_idx is None:
            num = 0
        else:
            num = last_idx + 1
        # increment twice for:
        # - add one to convert from zero-indexing to one-indexing
        # - row 1 is header
        row = num + 2
        return (prefix, num + 1), row

    def get_next_id(self, key):
        if key in self.ids:
            id_ = self.ids[key]
            (prefix, num), row = id_
            self.ids[key] = ((prefix, num + 1), row + 1)
            return id_
        sheet = self.get_sheet(key)

    def types_for_prefix(self, prefix):
        types = []
        for reg_prefix, type_ in self.registry.keys():
            if prefix == reg_prefix:
                types.append(type_)
        return types

    def get_map_list(self, prefix):
        if prefix in self.maps:
            return self.maps[prefix]
        else:
            key = (prefix, "maps")
            maps_folder = self.registry.get(key)
            if maps_folder is None:
                raise ValueError(f"expecting a registry entry {key}")
            maps = list_drive(
                self.sheets_client.drive.service, root=maps_folder, is_folder=False
            )
            self.maps[prefix] = maps
            return maps

    def _get_map(self, prefix, name):
        seq_files = self.get_map_list(prefix)
        seq_file = regex_key(seq_files, f"{name}\\.", check_duplicates=True)
        seq = read_sequence(
            self.sheets_client.drive.service.files()
            .get_media(fileId=seq_file["id"])
            .execute()
            .decode("utf8")
        )
        molecule_type = seq.annotations["molecule_type"]
        if molecule_type == "ds-DNA":
            seq = DsSeqRecord(seq, id=name, name=name, description=name)
        else:
            raise NotImplementedError(
                f"cannot import genbank molecule_type: '{molecule_type}'"
            )
        entry = {"_seq": seq}
        return entry

    def _get(self, df, prefix, type_, name):
        if name not in df.index:
            raise ValueError(f"cannot find {name}")
        entry = df.loc[name].to_dict()
        entry["_id"] = name
        entry["_type"] = type_
        entry["_prefix"] = prefix
        if type_ in TYPES_WITH_MAPS:
            entry.update(self._get_map(prefix, name))
        elif type_ == "parts":
            res = eval_exprs_by_priority(entry["Usage"], self.get)
            if res is None:
                seq = None
            else:
                seq = res.get("_seq")
            if seq is None:
                seq = part_entry_to_seq(entry)
            seq.id = seq.name = seq.description = name
            entry["_seq"] = seq
        elif "Sequence" in entry:
            entry["_seq"] = Seq(entry["Sequence"])
        return entry

    def get(self, name, types=("plasmids", "strains", "oligos", "parts"), seq=True):
        df, prefix, type_ = self.get_df_by_name(name, types=types)
        return self._get(df, prefix, type_, name)

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
