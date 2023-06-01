import re

import numpy as np
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


def nonempty_column_mask(worksheet, return_dataframe=False):
    df = worksheet.get_as_df(value_render=pygsheets.ValueRenderOption.FORMULA)
    has_datavalidation = columns_with_validation(
        worksheet.client.sheet.service,
        worksheet.spreadsheet.id,
        worksheet.spreadsheet._sheet_list,
    )
    mask = _nonempty_column_mask(df, has_datavalidation[worksheet.title])
    if return_dataframe:
        return mask, df
    else:
        return mask


def _formula_mask(df):
    return df.columns.str.startswith("=")


def _validation_mask(df, has_datavalidation):
    validation_mask = has_datavalidation
    validation_mask += [False] * (len(df.columns) - len(validation_mask))
    validation_mask = np.array(validation_mask)
    return validation_mask


def _nonempty_column_mask(df, has_datavalidation):
    formula_mask = _formula_mask(df)
    validation_mask = _validation_mask(df, has_datavalidation)
    mask = formula_mask | validation_mask
    return mask


###########


# def _get_next_empty_row(worksheet, skip_columns=0):
#     # TODO: we can probably remove these kwargs to get_as_df
#     # since pygsheets 2.0.4 fixed the bug with trailing empty cells
#     df = worksheet.get_as_df(has_header=False, empty_value=None)
#     df = df[1:]
#     has_datavalidation = columns_with_validation(
#         worksheet.client.sheet.service,
#         worksheet.spreadsheet.id,
#         worksheet.spreadsheet._sheet_list,
#     )
#     mask = has_datavalidation[worksheet.title]
#     mask += [False] * (len(df.columns) - len(mask))
#     nonempty = (
#         ~df.iloc[:, skip_columns:]
#         .iloc[:, ~np.array(mask)[skip_columns:]]
#         .isnull()
#         .all(axis=1)
#     )
#     last_idx = nonempty[nonempty].last_valid_index()
#     # convert to Python int because
#     # DataFrame.last_valid_index() returns np.int64, which is not JSON-serializable
#     last_idx = int(last_idx)
#     return last_idx, df


# def get_next_empty_row(worksheet, skip_columns=0):
#     last_idx, _ = _get_next_empty_row(worksheet, skip_columns=skip_columns)
#     if last_idx is None:
#         return 2
#     else:
#         # increment twice for:
#         # - add one to convert from zero-indexing to one-indexing
#         # - row 1 is header
#         return last_idx + 2


def get_next_empty_row(worksheet, skip_columns=0):
    mask, df = nonempty_column_mask(worksheet, return_dataframe=True)
    return _get_next_empty_row(df, mask, skip_columns=skip_columns)


def _get_next_empty_row(df, mask, skip_columns=0):
    masked_values = df.iloc[:, skip_columns:].iloc[:, ~mask[skip_columns:]]
    masked_values[masked_values == ""] = np.nan
    nonempty = ~masked_values.isnull().all(axis=1)
    last_idx = nonempty[nonempty].last_valid_index()
    if last_idx is not None:
        # convert to Python int because
        # DataFrame.last_valid_index() returns np.int64, which is not JSON-serializable
        last_idx = int(last_idx)
    return last_idx


def update_sheet_rows(service, updates, value_input_option="RAW"):
    if not updates:
        return
    body = {
        "value_input_option": value_input_option,
        "data": [
            {"range": datarange.range, "values": values}
            for datarange, values in updates
        ],
    }
    spreadsheet_ids = list(
        set(datarange.worksheet.spreadsheet.id for datarange, _ in updates)
    )
    if len(spreadsheet_ids) > 1:
        raise ValueError("all dataranges must refer to the same spreadsheet")
    request = (
        service.spreadsheets()
        .values()
        .batchUpdate(spreadsheetId=spreadsheet_ids[0], body=body)
    )
    request.execute()


def insert_sheet_rows(sheet, row, entries, ignore_missing_columns=False, columns=None):
    if columns is None:
        columns = sheet.get_row(1)
    seen_columns = set()
    values = []
    for entry in entries:
        seen_columns.update(entry.keys())
        values.append([entry.get(col, "") for col in columns])
    unused_columns = seen_columns - set(columns)
    if len(unused_columns) and not ignore_missing_columns:
        raise ValueError(
            f"attempting to insert rows with data for missing columns: {unused_columns}"
        )
    # we insert at row - 1 because this inserts below the given row number
    sheet.insert_rows(row - 1, number=len(values), values=values, inherit=True)


def clear_sheet(sheet, skiprows=1):
    first_row = sheet.get_row(
        skiprows + 1, value_render=pygsheets.ValueRenderOption.FORMULA
    )
    bounds = []
    for col, value in enumerate(first_row):
        if (
            not isinstance(value, str)
            or re.match(r"=HYPERLINK\(", value, re.IGNORECASE)
            or not value.startswith("=")
        ):
            bounds.append(((skiprows, col), (None, col + 1)))
    _clear_sheet_request(sheet, bounds)


def _clear_sheet_request(sheet, bounds, fields="userEnteredValue"):
    requests = []
    for start, end in bounds:
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
