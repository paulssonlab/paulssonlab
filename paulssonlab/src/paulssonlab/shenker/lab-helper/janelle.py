#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd
import click
from apiclient.discovery import build
from httplib2 import Http
from oauth2client import file, client, tools
import pendulum
from functools import wraps
from IPython import embed

CALENDAR_NAME = "Paulsson Lab Nikon and Accessories"
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


def get_calendar_service():
    SCOPES = "https://www.googleapis.com/auth/calendar.readonly"
    store = file.Storage("credentials.json")
    creds = store.get()
    if not creds or creds.invalid:
        flow = client.flow_from_clientsecrets("client_secret.json", SCOPES)
        creds = tools.run_flow(flow, store)
    service = build("calendar", "v3", http=creds.authorize(Http()))
    return service


def parse_time(time):
    # t = pendulum.parse(time.get('dateTime', time.get('date')), tz=time.get('timeZone'))
    t = pd.to_datetime(time.get("dateTime", time.get("date")))
    all_day = "date" in time and "dateTime" not in time
    return t, all_day


def iter_calendar_events(service, calendar_name, max_results=MAX_RESULTS):
    calendars = service.calendarList().list().execute()["items"]
    calendarId = next(cal["id"] for cal in calendars if cal["summary"] == calendar_name)
    for event in paginator(service.events().list)(
        calendarId=calendarId,
        singleEvents=True,
        orderBy="startTime",
        maxResults=max_results,
    ):
        yield event


def get_calendar_events(service, calendar_name, max_results=MAX_RESULTS):
    events = []
    for event in iter_calendar_events(service, calendar_name):
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
