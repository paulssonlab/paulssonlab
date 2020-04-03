import re
import io
from cytoolz import partial
from apiclient.http import MediaIoBaseUpload


def get_drive_id(url):
    return re.search(
        r"https://drive.google.com/drive/u/\d+/folders/(\S+)", url
    ).groups()[0]


def get_drive_query(service, query):
    response = service.files().list(q=query).execute()
    if "files" not in response or not response["files"]:
        raise ValueError(f"could not find file/folder '{p}'")
    files = response["files"]
    if len(files) > 1:
        raise ValueError(f"got multiple files/folders matching name '{p}'")
    return files[0]


def get_drive_by_name(service, name, folder=None):
    query = f"name = '{name}'"
    if folder is True:
        query += " and (mimeType = 'application/vnd.google-apps.folder')"
    elif folder is False:
        query += " and (mimeType != 'application/vnd.google-apps.folder')"
    return get_drive_query(service, query)["id"]


def get_drive_by_path(service, path, root=None, folder=None):
    path_components = path.split("/")
    for p in path_components[:-1]:
        query = f"(name = '{p}') and ('{root}' in parents) and (mimeType = 'application/vnd.google-apps.folder')"
        root = get_drive_query(service, query)["id"]
    query = f"(name = '{path_components[-1]}') and ('{root}' in parents)"
    if folder is True:
        query += " and (mimeType = 'application/vnd.google-apps.folder')"
    elif folder is False:
        query += " and (mimeType != 'application/vnd.google-apps.folder')"
    root = get_drive_query(service, query)["id"]
    return root


def filter_drive(files, keys):
    keys_to_id = {}
    for file in files:
        key = (file["name"], file["mimeType"] == "application/vnd.google-apps.folder")
        keys_to_id[key] = file["id"]
    filtered_files = {}
    for ref, key in keys.items():
        if key in keys_to_id:
            filtered_files[ref] = keys_to_id[key]
        else:
            raise ValueError(
                f"could not find {'folder' if key[1] else 'file'} '{key[0]}'"
            )
    return filtered_files


def list_drive(service, root=None, query=None):
    qs = []
    if root:
        qs.append(f"'{root}' in parents")
    if query:
        qs.append(query)
    q = " and ".join([f"({subq})" for subq in qs])
    response = service.files().list(q=q).execute()
    return _drive_list_to_dict(response)


def _drive_list_to_dict(response):
    files = {}
    if not "files" in response:
        return
    for file in response["files"]:
        if file["name"] in files:
            raise ValueError(f"got duplicate file name in drive list: {file['name']}")
        files[file["name"]] = file
    return files


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
