def paginator(req_func):
    @wraps(req_func)
    def f(*args, **kwargs):
        res = req_func(*args, **kwargs).execute()
        for item in res["items"]:
            yield item
        while "nextPageToken" in res and res["nextPageToken"]:
            res = req_func(
                *args, **{"pageToken": res["nextPageToken"], **kwargs}
            ).execute()
            for item in res["items"]:
                yield item

    return f


def get_creds():
    creds = None
    # The file token.pickle stores the user's access and refresh tokens, and is
    # created automatically when the authorization flow completes for the first
    # time.
    if os.path.exists("token.pickle"):
        with open("token.pickle", "rb") as token:
            creds = pickle.load(token)
    # If there are no (valid) credentials available, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file("credentials.json", SCOPES)
            creds = flow.run_local_server()
        # Save the credentials for the next run
        with open("token.pickle", "wb") as token:
            pickle.dump(creds, token)
    return creds


def get_calendar_service():
    service = build("calendar", "v3", credentials=get_creds())
    return service


def get_sheets_service():
    service = build("sheets", "v4", credentials=get_creds())
    return service


def get_calendar_id(service, calendar_name):
    calendars = service.calendarList().list().execute()["items"]
    calendar_id = next(
        cal["id"] for cal in calendars if cal["summary"] == calendar_name
    )
    return calendar_id


def iter_calendar_events(service, calendar_id, max_results=MAX_RESULTS, **kwargs):
    for event in paginator(service.events().list)(
        calendarId=calendar_id,
        singleEvents=True,
        orderBy="startTime",
        maxResults=max_results,
        **kwargs,
    ):
        yield event


def get_sheets_for_spreadsheet(sheets_service, spreadsheet_id):
    res = (
        sheets_service.spreadsheets()
        .get(spreadsheetId=spreadsheet_id, fields="sheets.properties")
        .execute()
    )
    sheets = {
        s["properties"]["title"]: s["properties"]["sheetId"] for s in res["sheets"]
    }
    return sheets
