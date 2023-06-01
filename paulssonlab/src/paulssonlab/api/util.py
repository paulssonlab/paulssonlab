import re

from tqdm.auto import tqdm

PROGRESS_BAR = tqdm
del tqdm


def _default_parse_html_table_row(column_names, tr):
    row = {}
    for name, td in zip(column_names, tr.find("td")):
        row[name] = td.text
    return row


def parse_html_table(table, row_parser=_default_parse_html_table_row):
    header = table.find("thead", first=True)
    rows = table.find("tr")
    column_names = [t.text for t in header.find("th")]
    rows = []
    trs = table.find("tbody tr")
    for tr in trs:
        row = row_parser(column_names, tr)
        rows.append(row)
    return rows


def base_url(url):
    return re.match("^(?:https?://)?(.*[^/]+)/?$", url).group(1).lower()


def regex_key(x, pattern, check_duplicates=True):
    pattern = re.compile(pattern)  # ensure it is compiled
    no_match = object()
    match = no_match
    for k in x:
        if pattern.match(k):
            returned_value = k
            if check_duplicates:
                if match == no_match:
                    match = returned_value
                else:
                    raise ValueError(
                        f"duplicate match found for pattern '{pattern.pattern}'"
                    )
            else:
                return returned_value
    if match == no_match:
        raise ValueError(f"no match found for pattern '{pattern.pattern}'")
    else:
        return match
