import re
import requests_html
import urllib
from paulssonlab.cloning.util import format_well_name
from paulssonlab.api.util import parse_html_table, PROGRESS_BAR


def _ensure_addgene_url(catalog_or_url):
    try:
        if "http" in catalog_or_url.lower():
            return catalog_or_url
    except:
        pass
    return f"https://www.addgene.org/{catalog_or_url}/"


def get_addgene_plasmid(catalog_or_url, session=None):
    url = _ensure_addgene_url(catalog_or_url)
    if not session:
        session = requests_html.HTMLSession()
    res = session.get(url)
    return _get_addgene_plasmid(res.html)


def _get_addgene_plasmid(html):
    material = {}
    # name ("pAJM.711")
    material["name"] = html.find("div.panel-heading span.material-name")[0].text
    # add-to-cart table
    table = parse_html_table(html.find("table.add-to-cart-table", first=True))
    for row in table:
        del row[""]
    if len(table) != 1:
        raise ValueError("expecting a single row in Addgene add-to-cart table")
    # item ("Plasmid")
    material["item"] = table[0]["Item"]
    # catalog ("108512")
    catalog_number = int(table[0]["Catalog #"])
    material["catalog"] = catalog_number
    material["url"] = f"https://www.addgene.org/{catalog_number}/"
    # price
    material["price"] = table[0]["Price (USD)"]
    # purpose, depositing lab (text -> href), publication (text -> href)
    lis = html.find("ul#plasmid-description-list > li")
    for li in lis:
        label = li.find("div.field-label", first=True)
        if not label:
            continue
        key = label.text
        content = li.find("div.field-content", first=True)
        if not content:
            continue
        link = content.find("a", first=True)
        if link:
            value = {link.text: urllib.parse.urljoin(link.url, link.attrs["href"])}
        else:
            value = content.text
        material[key.lower()] = value
    # MAPPING:
    # bacterial resistance(s)
    # growth temperature
    # growth strain
    # copy number
    # gene insert name
    # gene/insert species
    # cloning method
    # supp docs (text -> href)
    # licenses (text -> href)
    # depositor comments
    sections = html.find("div#detail-sections section")
    for section in sections:
        title = section.find("h2 > span.title", first=True)
        key = title.text.lower()
        content = section.find("h2 + div", first=True)
        if content:
            material[key] = content.text
        else:
            ul = section.find("h2 + ul", first=True)
            if not ul:
                continue
            # to get immediate children of ul, we need to use PyQuery directly in this hacky way
            for li in ul.pq.children("ul > li"):
                # re-wrap lxml.html.HtmlElement in html_requests.Element
                li = requests_html.Element(
                    element=li, url=ul.url, default_encoding=ul.default_encoding
                )
                label = li.find(".field-label", first=True)
                key = label.text.lower()
                doc_list = li.find("ul.addgene-document-list", first=True)
                if doc_list:
                    doc_links = doc_list.find("li > a")
                    value = {
                        link.text: urllib.parse.urljoin(link.url, link.attrs["href"])
                        for link in doc_links
                    }
                else:
                    value = label.element.tail.strip()
                material[key] = value
    # how to cite (materials_and_methods)
    # how to cite (references)
    how_to_cite = {}
    how_to_cite_lis = html.find("section#how-to-cite li")
    for li in how_to_cite_lis:
        key = li.find("strong", first=True).text.lower().replace(" & ", "_and_")
        value = li.find("small", first=True).text
        how_to_cite[key] = value
    material["how_to_cite"] = how_to_cite
    return material


def get_addgene_sequence_urls(url, session=None):
    if not session:
        session = requests_html.HTMLSession()
    res = session.get(urllib.parse.urljoin(url, "sequences"))
    return _get_addgene_sequence_urls(res.html)


def _get_addgene_sequence_urls(html):
    seqs = {}
    for key in [
        "addgene_full",
        "depositor_full",
        "addgene_partial",
        "depositor_partial",
    ]:
        links = html.find(f"section#{key.replace('_', '-')} a.genbank-file-download")
        if links:
            seq_urls = [link.attrs["href"] for link in links]
            seqs[key] = seq_urls
    return seqs


def _parse_addgene_well(s):
    m = re.match(r"(?:Plate\s+(\d+) / )?([A-H]) / (\d+)", s.strip())
    return m.groups()


def _parse_addgene_kit_row(column_names, tr):
    url = None
    tds = tr.find("td")
    for name, td in zip(column_names, tds):
        link = td.find("a", first=True)
        if link is not None and not link.attrs["href"].startswith("#"):
            url = urllib.parse.urljoin(link.base_url, link.attrs["href"])
    row = {name.lower(): td.text for name, td in zip(column_names, tds)}
    if url is not None:
        row["url"] = url
        row["catalog"] = int(
            re.match(r"https?://www.addgene.org/(\d+)/?", url).group(1)
        )
    if "well" in row:
        row["well"] = format_well_name(*_parse_addgene_well(row["well"]))
    return row


def get_addgene(
    catalog_or_url,
    include_sequences=True,
    include_details=False,
    progress_bar=PROGRESS_BAR,
    session=None,
):
    url = _ensure_addgene_url(catalog_or_url)
    if not session:
        session = requests_html.HTMLSession()
    res = session.get(url)
    return _get_addgene(
        res.html,
        include_sequences=include_sequences,
        include_details=include_details,
        progress_bar=progress_bar,
        session=session,
    )


def _get_addgene(
    html,
    include_sequences=True,
    include_details=False,
    progress_bar=PROGRESS_BAR,
    session=None,
):
    data = {}
    catalog = html.find("small#catalog-number", first=True)
    if catalog is None:
        catalog = html.find("span.material-name ~ small", first=True)
    item, catalog_number = re.search(
        r"\(([\w ]+) #\s*(\d+)\s*\)", catalog.text
    ).groups()
    catalog_number = int(catalog_number)
    data["item"] = item
    data["catalog"] = catalog_number
    url = f"http://www.addgene.org/{catalog_number}/"
    data["url"] = url
    if item == "Kit":
        kit = _get_addgene_kit(
            html,
            include_sequences=include_sequences,
            include_details=include_details,
            progress_bar=progress_bar,
        )
        return {**data, **kit}
    elif item == "Plasmid" or item == "Bacterial strain":
        plasmid = _get_addgene_plasmid(html)
        if include_sequences:
            sequence_urls = get_addgene_sequence_urls(url, session=session)
            plasmid["sequence_urls"] = sequence_urls
        return plasmid
    else:
        raise ValueError(f"unknown Addgene item type: {item}")


def _get_addgene_kit(
    html,
    include_sequences=True,
    include_details=False,
    progress_bar=PROGRESS_BAR,
    session=None,
):
    table = html.find("table.kit-inventory-table", first=True)
    if table is None:
        table = html.find("div#kit-contents > table", first=True)
    if table is None:
        raise ValueError("could not find Addgene kit table")
    wells = parse_html_table(table, row_parser=_parse_addgene_kit_row)
    if progress_bar is not None and (include_sequences or include_details):
        wells_iter = progress_bar(wells)
    else:
        wells_iter = wells
    for well in wells_iter:
        if include_sequences:
            sequence_urls = get_addgene_sequence_urls(well["url"], session=session)
            well["sequence_urls"] = sequence_urls
        if include_details:
            plasmid = get_addgene_plasmid(well["url"], session=session)
            well.update(plasmid)
    data = {"wells": wells}
    return data
