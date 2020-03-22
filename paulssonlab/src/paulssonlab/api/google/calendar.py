MAX_RESULTS = 2500


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
