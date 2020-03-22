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
