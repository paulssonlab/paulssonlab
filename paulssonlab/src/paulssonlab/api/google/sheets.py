import pygsheets

SHEETS_MIMETYPE = "application/vnd.google-apps.spreadsheet"


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


def columns_with_validation(service, spreadsheet_id, worksheets=None):
    if worksheets is not None:
        ranges = [
            f"'{w.title}'!2:2" for w in worksheets
        ]  # get the row after the header
    else:
        ranges = "2:2"
    params = {
        "spreadsheetId": spreadsheet_id,
        "ranges": ranges,
        "fields": "sheets(data/rowData/values/dataValidation,properties(sheetId,title))",
    }
    response = service.spreadsheets().get(**params).execute()
    has_datavalidation = {}
    for sheet in response["sheets"]:
        sheet_name = sheet["properties"]["title"]
        if "data" in sheet and "rowData" in sheet["data"][0]:
            datavalidation = sheet["data"][0]["rowData"][0]["values"]
            has_datavalidation[sheet_name] = [
                "dataValidation" in col for col in datavalidation
            ]
        else:
            has_datavalidation[sheet_name] = []
    return has_datavalidation


def insert_sheet_rows(sheet, row, entries):
    columns = sheet.get_row(1)
    values = [[entry.get(col) for col in columns] for entry in entries]
    # we insert at row - 1 because this inserts below the given row number
    sheet.insert_rows(row - 1, number=len(values), values=values)


def clear_sheet(sheet, skiprows=1):
    first_row = sheet.get_row(
        skiprows + 1, value_render=pygsheets.ValueRenderOption.FORMULA
    )
    bounds = []
    for col, value in enumerate(first_row):
        if not isinstance(value, str) or not value.startswith("="):
            bounds.append(((skiprows, col), (None, col + 1)))
    _clear_sheet_request(sheet, bounds)


def _clear_sheet_request(sheet, bounds, fields="userEnteredValue"):
    requests = []
    for (start, end) in bounds:
        range_ = {
            "sheetId": sheet.id,
            "startRowIndex": start[0],
            "startColumnIndex": start[1],
            "endRowIndex": end[0],
            "endColumnIndex": end[1],
        }
        request = {"updateCells": {"range": range_, "fields": fields}}
        requests.append(request)
    sheet.client.sheet.batch_update(sheet.spreadsheet.id, requests)
