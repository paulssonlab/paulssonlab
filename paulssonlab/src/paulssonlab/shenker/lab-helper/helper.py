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
from datetime import datetime
from functools import wraps
from itertools import zip_longest
from collections import defaultdict
import re
from IPython import embed

# TODO: individual meeting duration (45min)
# TODO: move config to .toml
SCOPES = [
    "https://www.googleapis.com/auth/calendar.events",
    "https://www.googleapis.com/auth/calendar.readonly",
    "https://www.googleapis.com/auth/spreadsheets",
]
GENERAL_CALENDAR = "Paulsson Lab General Calendar"
MICROSCOPE_CALENDAR = "Paulsson Lab Nikon and Accessories"
GROUP_MEETINGS_SPREADSHEET = "1Wm3U7YZHewkRthGeDZJ-ahqAY5FW8w91QkAoqsid1aQ"
INDIVIDUAL_MEETINGS_SPREADSHEET = "18ljx4U1RednlX555YSFbQWr0ziz2FmwadKdut358qk0"
GROUP_MEETING_DURATION = 2
INDIVIDUAL_MEETING_DURATION = 45 / 60
GOOGLE_USER = "paulssonlab@gmail.com"
MAX_RESULTS = 2500
DEFAULT_TIMEZONE = "America/New_York"


def none_if_exception(func):
    try:
        return func()
    except:
        return None


def parse_date_and_time(date, time, tz=pendulum.timezone("utc")):
    # date = pendulum.parse(date, exact=True)
    # time = pendulum.parse(time, exact=True)
    # FROM: https://github.com/sdispater/pendulum/issues/156
    # return pendulum.instance(datetime.combine(date, time), tz=tz)
    return pendulum.parse(date + " " + time, tz=tz, strict=False)


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


def parse_time(time, all_day_tz=DEFAULT_TIMEZONE):
    # t = pendulum.parse(time.get('dateTime', time.get('date')), tz=time.get('timeZone'))
    # t = pd.to_datetime(time.get('dateTime', time.get('date')))
    all_day = "date" in time and "dateTime" not in time
    if all_day:
        t = pd.to_datetime(time["date"])
        t = t.tz_localize(all_day_tz).tz_convert("UTC")
    else:
        t = pd.to_datetime(time.get("dateTime", time.get("date")))
        if not t.tz:
            t = t.tz_localize(time["timeZone"])
        t = t.tz_convert("UTC")
    return t, all_day


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


def get_calendar_events(
    service, calendar_id, max_results=MAX_RESULTS, all_day_tz=DEFAULT_TIMEZONE, **kwargs
):
    events = []
    for event in iter_calendar_events(service, calendar_id, **kwargs):
        start, all_day = parse_time(event["start"], all_day_tz=all_day_tz)
        end, _ = parse_time(event["end"], all_day_tz=all_day_tz)
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
    service, calendar_id, max_results=MAX_RESULTS, **kwargs
):
    now = pendulum.now().isoformat()
    events = []
    for event in iter_calendar_events(
        service,
        calendar_id,
        timeMin=now,
        sharedExtendedProperty="sync=lab-helper",
        max_results=max_results,
        **kwargs,
    ):
        events.append(event)
    return events


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


def _get_sheets_group_meetings(spreadsheet_id=GROUP_MEETINGS_SPREADSHEET):
    sheets_service = get_sheets_service()
    res = (
        sheets_service.spreadsheets()
        .values()
        .get(spreadsheetId=spreadsheet_id, range="A:Z")
        .execute()
    )
    columns = res["values"][0]
    meetings = [
        {col: d.strip() for col, d in zip_longest(columns, row, fillvalue="")}
        for row in res["values"][1:]
    ]
    return meetings


def get_sheets_group_meetings(spreadsheet_id=GROUP_MEETINGS_SPREADSHEET):
    meetings = []
    for m in _get_sheets_group_meetings(spreadsheet_id=spreadsheet_id):
        summary, description = _format_group_meeting(m)
        del m["Subject"]
        del m["Cook"]
        del m["Notes"]
        m["summary"] = summary
        m["description"] = description
        meetings.append(m)
    return meetings


def get_sheets_individual_meetings(spreadsheet_id=INDIVIDUAL_MEETINGS_SPREADSHEET):
    sheets_service = get_sheets_service()
    sheets = get_sheets_for_spreadsheet(sheets_service, spreadsheet_id)
    meetings = []
    for sheet_title in sheets:
        # TODO: get sheet using sheetId not sheet_title
        # (there might be escaping issues with sheet_title)
        # SEE: https://stackoverflow.com/questions/43667019/looking-up-google-sheets-values-using-sheet-id-instead-of-sheet-title-google-dr/43837549
        # SEE: https://stackoverflow.com/questions/43916341/php-google-sheets-api-v4-not-working-with-sheet-names-containing-single-quotes
        res = (
            sheets_service.spreadsheets()
            .values()
            .get(
                spreadsheetId=spreadsheet_id,
                range="'{}'".format(sheet_title),
                majorDimension="COLUMNS",
            )
            .execute()
        )
        time_col = res["values"][0]
        idx = next((i for i, x in enumerate(time_col) if x), None)
        times = time_col[idx:]
        dates = [col[idx - 1] for col in res["values"][1:]]
        for date, col in zip(dates, res["values"][1:]):
            for time, name in zip(times, col[idx:]):
                if not name or "overflow" in name:
                    continue
                meeting = {
                    "Date": date,
                    "Time": time,
                    "Room": "Johan's office",
                    "summary": "Johan: {}".format(name),
                    "description": "",
                }
                meetings.append(meeting)
    return meetings


# @click.group()
# @click.pass_context
# def cli(ctx):
#     service = get_calendar_service()
#     cfg = {'service': service}
#     ctx.obj = cfg

# @cli.command()
# @click.pass_obj
# def test(cfg):
#     service = cfg['service']
#     embed()


@click.group()
def cli():
    pass


def _format_group_meeting(m):
    # only create events in the future
    # TODO: tbd, canceled, czar, empty fields, notes
    if not m["Speaker"] and not m["Subject"]:
        return None, None
    if "canceled" in m["Subject"].casefold():
        return "CANCELED: group meeting", m["Notes"]
    if "(czar)" in m["Speaker"].casefold() or "subgroup" in m["Subject"].casefold():
        # subgroup meeting
        speaker = re.sub(r"\s*\(?czar\)?\s*", "", m["Speaker"], re.IGNORECASE)
        summary = "{} (czar: {})".format(m["Subject"], speaker)
    else:
        if m["Speaker"]:
            summary = m["Speaker"]
            if m["Subject"]:
                summary += ": {}".format(m["Subject"])
        else:
            summary = m["Subject"]
    descriptions = []
    if m["Cook"]:
        descriptions.append("Cook: {}".format(m["Cook"]))
    descriptions.append(m["Notes"])
    description = "\n\n".join(descriptions)
    return summary, description


@cli.command()
@click.option("--group", default=True)
@click.option("--individual", default=False)
def sync_meetings(group, individual):
    calendar_service = get_calendar_service()
    calendar_id = get_calendar_id(calendar_service, GENERAL_CALENDAR)
    sheets_meetings = []
    if group:
        for m in get_sheets_group_meetings():
            m["duration"] = GROUP_MEETING_DURATION
            sheets_meetings.append(m)
    if individual:
        for m in get_sheets_individual_meetings():
            m["duration"] = INDIVIDUAL_MEETING_DURATION
            sheets_meetings.append(m)
    old_events = get_syncable_calendar_events(calendar_service, calendar_id)
    old_events_by_datetime = {e["start"]["dateTime"]: e for e in old_events}
    event_ids_to_update = set()
    new_events = []
    now = pendulum.now()
    for m in sheets_meetings:
        if m["summary"] is None and m["description"] is None:
            continue
        start_time = parse_date_and_time(
            m["Date"], m["Time"], tz=pendulum.timezone(TIMEZONE)
        )
        if start_time < now:
            # only add/update events in the future
            continue
        end_time = start_time.add(hours=m["duration"])
        event = {
            "summary": m["summary"],
            "location": m["Room"],
            "description": m["description"],
            "start": {"dateTime": start_time.isoformat()},
            "end": {"dateTime": end_time.isoformat()},
            "extendedProperties": {"shared": {"sync": "lab-helper"}},
        }
        new_events.append(event)
    for new_event in new_events:
        event_id_to_update = dict.get(
            old_events_by_datetime.get(new_event["start"]["dateTime"], None) or {}, "id"
        )
        if event_id_to_update:
            event_ids_to_update.add(event_id_to_update)
            calendar_service.events().update(
                calendarId=calendar_id, eventId=event_id_to_update, body=new_event
            ).execute()
        else:
            calendar_service.events().insert(
                calendarId=calendar_id, body=new_event
            ).execute()
    for event in old_events:
        if event["id"] not in event_ids_to_update:
            calendar_service.events().delete(
                calendarId=calendar_id, eventId=event["id"]
            ).execute()


if __name__ == "__main__":
    cli()
