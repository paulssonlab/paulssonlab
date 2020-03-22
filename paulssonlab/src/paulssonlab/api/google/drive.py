import re


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
