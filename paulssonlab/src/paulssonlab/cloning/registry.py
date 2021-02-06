import re
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
)
from paulssonlab.api import read_sequence, regex_key
from paulssonlab.cloning.sequence import DsSeqRecord, anneal, pcr
from paulssonlab.cloning.commands import expr_parser, command_parser

DEFAULT_LOOKUP_TYPES = ["oligos", "plasmids", "strains", "parts"]
FOLDER_TYPES = ["maps", "sequencing"]
TYPES_WITH_SEQUENCING = ["plasmids", "strains"]
TYPES_WITH_MAPS = ["plasmids"]
TYPES_WITHOUT_IDS = ["parts"]
COLLECTION_REGEX = r"^([^_]+)_(\w+)$"
ENABLE_AUTOMATION_FILENAME = "ENABLE_AUTOMATION.txt"
EXPR_PRIORITY = ["digest", "pcr"]


def _name_mapper_for_prefix(old_prefix, new_prefix):
    mapper = lambda name, is_folder: re.sub(
        rf"(\w*){re.escape(old_prefix)}(.*)", rf"\1{new_prefix}\2", name
    )
    return mapper


class Registry(object):
    def __init__(self, sheets_client, registry_folder):
        self.sheets = {}
        self.dfs = {}
        self.maps = {}
        self.sheets_client = sheets_client
        self.registry_folder = registry_folder
        self.refresh()

    def clear_cache(self):
        self.sheets = {}
        self.dfs = {}
        self.maps = {}

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

    def get_df_by_id(self, id_):
        df = self.dfs.get(id_)
        if df is None:
            # include_tailing_empty=True is necessary to prevent an error if the last column of data is empty
            # TODO: this may no longer be true in pygsheets 2.0.4
            df = self.get_sheet_by_id(id_).get_as_df(
                index_column=1, include_tailing_empty=True
            )
            self.dfs[id_] = df
        return df

    def get_df(self, name, types=("plasmids", "strains", "oligos", "parts")):
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
        seq_file = regex_key(seq_files, f"{name}\\.", check_duplicates=True)["id"]
        seq = read_sequence(
            self.sheets_client.drive.service.files()
            .get_media(fileId=seq_file)
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
        return seq

    def eval_exprs(self, s):
        if not s.strip():
            return None
        ast = expr_parser.parse(s)
        expr = None
        for type_ in EXPR_PRIORITY:
            priority_expr = [e for e in ast if e["_type"] == type_]
            if len(priority_expr):
                expr = priority_expr[0]
                break
        if expr is None and len(ast):
            expr = ast[0]
        return self.eval_expr(expr)

    def eval_expr(self, expr):
        if expr is None:
            return None
        type_ = expr["_type"]
        if type_ == "pcr":
            return pcr(
                self.eval_expr(expr["template"]),
                self.eval_expr(expr["primer1"]),
                self.eval_expr(expr["primer2"]),
            )
        elif type_ == "digest":
            enzyme_name = expr["enzyme"]["name"]
            if not hasattr(Bio.Restriction, enzyme_name):
                raise ValueError(f"unknown enzyme '{enzyme_name}'")
            enzyme = getattr(Bio.Restriction, enzyme_name)
            return re_digest_part(self.eval_expr(expr["input"]), enzyme)
        elif type_ == "anneal":
            return anneal(
                self.eval_expr(expr["strand1"]), self.eval_expr(expr["strand2"])
            )
        elif type_ == "name":
            entry = self.get(expr["name"])
            if entry is None or "_seq" not in entry:
                return None
            else:
                return entry["_seq"]
        else:
            return NotImplementedError

    def command(self, s):
        pass

    def get(self, name, types=("plasmids", "strains", "oligos", "parts"), seq=True):
        df, prefix, type_ = self.get_df(name, types=types)
        if name not in df.index:
            raise ValueError(f"cannot find {name}")
        entry = df.loc[name].to_dict()
        entry["_id"] = name
        entry["_type"] = type_
        entry["_prefix"] = prefix
        if type_ in TYPES_WITH_MAPS:
            entry["_seq"] = self._get_map(prefix, name)
        elif type_ == "parts":
            seq = self.eval_exprs(entry["Usage"])
            if seq is None:
                seq = part_entry_to_seq(entry)
            entry["_seq"] = seq
        elif "Sequence" in entry:
            entry["_seq"] = Seq(entry["Sequence"])
        return entry
