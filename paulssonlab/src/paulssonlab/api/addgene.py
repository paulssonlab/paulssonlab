import requests_html
import urllib
from paulssonlab.cloning.util import format_well_name


def _parse_addgene_well(s):
    m = re.match(r"Plate (\d+) / ([A-H]) / (\d+)", s)
    return m.groups()


def _parse_table_row(column_names, tr):
    url = None
    for name, td in zip(column_names, tr.find("td")):
        link = td.find("a", first=True)
        if link is not None and not link.attrs["href"].startswith("#"):
            url = urllib.parse.urljoin(link.base_url, link.attrs["href"])
    row = {name: td.text for name, td in zip(column_names, tr.find("td"))}
    if url is not None:
        row["url"] = url
    if "Well" in row:
        row["Well"] = format_well_name(*_parse_addgene_well(row["Well"]))
    return d


def addgene_supplemental_urls(url, session=None):
    if not session:
        session = requests_html.HTMLSession()
    res = session.get(url)
    return _addgene_supplemental_urls(res.html)


def _addgene_supplemental_urls(html):
    urls = html.find(
        "div.field-label:contains('Supplemental') + ul.addgene-document-list"
    )[0].links
    return urls


def addgene_sequences(url, session=None):
    if not session:
        session = requests_html.HTMLSession()
    res = session.get(urllib.parse.urljoin(url, "sequences"))
    return _addgene_sequences(res.html)


def _addgene_sequences(html):
    seqs = {}
    for key in [
        "addgene_full",
        "depositor_full",
        "addgene_partial",
        "depositor_partial",
    ]:
        links = html.find(f"section#{key.replace('_', '-')} a.genbank-file-download")
        seq_urls = [link.attrs["href"] for link in links]
        seqs[key] = seq_urls
    return seqs


def parse_addgene_kit_table(url, include_sequences=True, include_supplemental=True):
    session = requests_html.HTMLSession()
    res = session.get(url)
    table = res.html.find("table.kit-inventory-table")[0]
    header = table.find("thead", first=True)
    rows = table.find("tr")
    column_names = [t.text for t in header.find("th")]
    wells = []
    for row in table.find("tbody tr"):
        well = _parse_table_row(column_names, row)
        if include_sequences:
            sequence_urls = addgene_supplemental_urls(row["url"], session=session)
            well["sequence_urls"] = sequence_urls
        if include_supplemental:
            supp_urls = addgene_supplemental_links(row["url"], session=session)
            well["supplemental_urls"] = supp_urls
        wells.append(well)
    return wells
