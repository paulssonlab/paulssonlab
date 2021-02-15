import re
import io
from cytoolz import partial
from apiclient.http import MediaIoBaseUpload
from google.api_core.datetime_helpers import from_rfc3339

FOLDER_MIMETYPE = "application/vnd.google-apps.folder"


def get_drive_id(url):
    return re.search(
        r"https://drive.google.com/drive/u/\d+/folders/(\S+)", url
    ).groups()[0]


def get_drive_query(service, query):
    response = service.files().list(q=query).execute()
    if "files" not in response or not response["files"]:
        raise ValueError(f"could not find file/folder with query '{query}'")
    files = response["files"]
    if len(files) > 1:
        raise ValueError(f"got multiple files/folders matching query '{query}'")
    return files[0]


def get_drive_by_path(service, path, root=None, is_folder=None):
    if isinstance(path, str):
        path = [path]
    for p in path[:-1]:
        query = f"(name = '{p}') and (mimeType = '{FOLDER_MIMETYPE}')"
        if root is not None:
            query += f" and ('{root}' in parents)"
        root = get_drive_query(service, query)["id"]
    query = f"(name = '{path[-1]}')"
    if root is not None:
        query += f" and ('{root}' in parents)"
    if is_folder is True:
        query += f" and (mimeType = '{FOLDER_MIMETYPE}')"
    elif is_folder is False:
        query += f" and (mimeType != '{FOLDER_MIMETYPE}')"
    root = get_drive_query(service, query)["id"]
    return root


def list_drive(
    service, root=None, query=None, is_folder=None, page_size=1000, fields=[]
):
    qs = []
    if root:
        qs.append(f"'{root}' in parents")
    if is_folder is True:
        qs.append(f"mimeType = '{FOLDER_MIMETYPE}'")
    elif is_folder is False:
        qs.append(f"mimeType != '{FOLDER_MIMETYPE}'")
    if query:
        qs.append(query)
    q = " and ".join([f"({subq})" for subq in qs])
    files_service = service.files()
    fields = ["kind", "id", "name", "mimeType", *fields]
    fields = "files({})".format(",".join(fields))
    req = files_service.list(q=q, pageSize=page_size, fields=fields)
    files = []
    while req is not None:
        doc = req.execute()
        files.extend(doc.get("files", []))
        req = files_service.list_next(req, doc)
    name_to_file = {}
    for file in files:
        if file["name"] in name_to_file:
            raise ValueError(f"got duplicate file name in drive list: {file['name']}")
        name_to_file[file["name"]] = file
    for file in name_to_file.values():
        if "modifiedTime" in file:
            file["modifiedTime"] = from_rfc3339(file["modifiedTime"])
    return name_to_file


def ensure_drive_folder(file, is_folder=True):
    # could also use file["kind"] == "drive#folder"
    if (file["mimeType"] == FOLDER_MIMETYPE) != is_folder:
        raise ValueError(
            f"expecting {'folder' if is_folder else 'file'} for file '{file['name']}'"
        )
    return file["id"]


def make_drive_folder(service, name, parent):
    new_folder_metadata = {
        "name": name,
        "mimeType": FOLDER_MIMETYPE,
        "parents": [parent],
    }
    new_folder = service.files().create(body=new_folder_metadata, fields="id").execute()
    return new_folder["id"]


def copy_drive_file(service, source, name, parent):
    new_body = {"name": name, "parents": [parent]}
    new_file = service.files().copy(fileId=source, body=new_body).execute()
    return new_file["id"]


def copy_drive_folder(
    service,
    source_folder,
    dest_folder,
    folders_only=False,
    transform_names=None,
    recursive=True,
    overwrite=True,
):
    source_files = list_drive(service, source_folder)
    dest_files = list_drive(service, dest_folder)
    for source_file in source_files.values():
        if transform_names:
            new_name = transform_names(
                source_file["name"], source_file["mimeType"] == FOLDER_MIMETYPE
            )
        else:
            new_name = source_file["name"]
        if source_file["mimeType"] != FOLDER_MIMETYPE:
            # file
            if folders_only:
                continue
            if new_name in dest_files:
                if overwrite:
                    existing_file = dest_files[new_name]["id"]
                    service.files().delete(fileId=existing_file).execute()
                else:
                    raise ValueError(
                        f"attempting to overwrite existing file: '{new_name}' (id: {dest_files[new_name]['id']})"
                    )
            copy_drive_file(service, source_file["id"], new_name, dest_folder)
        else:
            # folder
            folder = None
            if new_name in dest_files:
                if overwrite:
                    folder = dest_files[new_name]["id"]
                    service.files().delete(fileId=existing_folder).execute()
                else:
                    if dest_files[new_name]["mimeType"] != FOLDER_MIMETYPE:
                        raise ValueError(
                            f"expecting : '{new_name}' (id: {dest_files[new_name]['id']})"
                        )
                    # if there's an existing folder with the right name, just use it
                folder = dest_files[new_name]["id"]
            if folder is None:
                folder = make_drive_folder(service, new_name, dest_folder)
            if recursive:
                copy_drive_folder(
                    service,
                    source_file["id"],
                    folder,
                    folders_only=folders_only,
                    transform_names=transform_names,
                    recursive=recursive,
                    overwrite=overwrite,
                )


def get_drive_modified_time(service, file_id):
    res = service.files().get(fileId=file_id, fields="modifiedTime").execute()
    modified_time = res.get("modifiedTime")
    if modified_time:
        return from_rfc3339(modified_time)


def upload_drive(
    service,
    content,
    name=None,
    mimetype="application/octet-stream",
    file_id=None,
    parent=None,
):
    if not hasattr(content, "seek"):
        content = io.BytesIO(content.encode())
    media = MediaIoBaseUpload(content, mimetype=mimetype, resumable=True)
    body = {"mimeType": mimetype}
    if name is not None:
        body["name"] = name
    if file_id is not None:
        method = partial(service.files().update, fileId=file_id)
    else:
        if name is None:
            raise ValueError("either file_id or name must be given")
        if parent is not None:
            body["parents"] = [parent]
        method = service.files().create
    response = method(body=body, media_body=media, fields="id").execute()
    return response.get("id")
