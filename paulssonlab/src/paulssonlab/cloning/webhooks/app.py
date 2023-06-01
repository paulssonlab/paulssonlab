import os
from datetime import datetime

import flask
import google.auth

# from googleapiclient.discovery import build
import pygsheets

# SCOPES = ["https://www.googleapis.com/auth/spreadsheets"]
SPREADSHEET_ID = "1Ab1h1_WzcJ4MMGPwl6qNwbXWMCQGG3ML4bpnGXvN1yU"

app = flask.Flask(__name__)

# initialize GCP credentials on cold start
credentials, project_id = google.auth.default()
gc = pygsheets.client.Client(credentials)


@app.route("/", methods=["GET", "POST"])
def hello_world():
    # return "Creds: {}\nProject ID: {}".format(type(credentials), project_id)
    # target = os.environ.get("TARGET", "World")
    # service = build("sheets", "v4")
    # sheet = service.spreadsheets()
    # result = sheet.values().get(spreadsheetId=SPREADSHEET_ID, range="A1:A").execute()
    # values = result.get("values", [])
    # gc = pygsheets.authorize()
    spreadsheet = gc.open_by_key(SPREADSHEET_ID)
    worksheet = spreadsheet[0]
    # old_val = worksheet.get_value("A1")
    val = str(flask.request.headers)
    # val = "Now: {}".format(datetime.now().isoformat())
    worksheet.update_value("A1", val)
    # return "Hello 3 {}!\n".format(old_val)
    return "Hello"


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))
