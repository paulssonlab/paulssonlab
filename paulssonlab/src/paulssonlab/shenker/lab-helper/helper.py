#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd
import click
import pickle
import os.path
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request
import pendulum
from functools import wraps
from itertools import zip_longest
from IPython import embed

SCOPES = [
    "https://www.googleapis.com/auth/calendar.events",
    "https://www.googleapis.com/auth/calendar.readonly",
    "https://www.googleapis.com/auth/spreadsheets",
]
GENERAL_CALENDAR = "Paulsson Lab General Calendar"
MICROSCOPE_CALENDAR = "Paulsson Lab Nikon and Accessories"
GROUP_MEETINGS_SHEET = "1Wm3U7YZHewkRthGeDZJ-ahqAY5FW8w91QkAoqsid1aQ"
INDIVIDUAL_MEETINGS_SHEET = "18ljx4U1RednlX555YSFbQWr0ziz2FmwadKdut358qk0"
GOOGLE_USER = "paulssonlab@gmail.com"
MAX_RESULTS = 2500


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


def none_if_exception(func):
    try:
        return func()
    except:
        return None


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


def parse_time(time):
    # t = pendulum.parse(time.get('dateTime', time.get('date')), tz=time.get('timeZone'))
    t = pd.to_datetime(time.get("dateTime", time.get("date")))
    all_day = "date" in time and "dateTime" not in time
    return t, all_day


def iter_calendar_events(service, calendar_name, max_results=MAX_RESULTS, **kwargs):
    calendars = service.calendarList().list().execute()["items"]
    calendarId = next(cal["id"] for cal in calendars if cal["summary"] == calendar_name)
    for event in paginator(service.events().list)(
        calendarId=calendarId,
        singleEvents=True,
        orderBy="startTime",
        maxResults=max_results,
        **kwargs,
    ):
        yield event


def get_calendar_events(service, calendar_name, max_results=MAX_RESULTS, **kwargs):
    events = []
    for event in iter_calendar_events(service, calendar_name, **kwargs):
        start, all_day = parse_time(event["start"])
        end, _ = parse_time(event["end"])
        created = none_if_exception(lambda: pd.to_datetime(event["created"]))
        updated = none_if_exception(lambda: pd.to_datetime(event["updated"]))
        events.append(
            {
                "summary": event["summary"],
                "start": start,
                "end": end,
                "all_day": all_day,
                "creator": event["creator"].get(
                    "displayName", event["creator"]["email"]
                ),
                "created": created,
                "updated": updated,
            }
        )
    return pd.DataFrame(events)


def get_syncable_calendar_events(
    service, calendar_name, max_results=MAX_RESULTS, **kwargs
):
    now = pendulum.now().isoformat()
    events = []
    for event in iter_calendar_events(
        service,
        calendar_name,
        timeMin=now,
        sharedExtendedProperty="sync=janelle",
        max_results=max_results,
        **kwargs,
    ):
        events.append(event)
    return events


def get_group_meeting_list():
    sheets_service = get_sheets_service().spreadsheets()
    res = (
        sheets_service.values()
        .get(spreadsheetId=GROUP_MEETINGS_SHEET, range="A:Z")
        .execute()
    )
    columns = res["values"][0]
    meetings = [
        {col: d for col, d in zip_longest(columns, row, fillvalue="")}
        for row in res["values"][1:]
    ]
    return meetings


@click.group()
@click.pass_context
def cli(ctx):
    service = get_calendar_service()
    cfg = {"service": service}
    ctx.obj = cfg


@cli.command()
@click.pass_obj
def test(cfg):
    service = cfg["service"]

    embed()


if __name__ == "__main__":
    cli()
