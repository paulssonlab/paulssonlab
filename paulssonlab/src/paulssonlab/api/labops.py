import pandas as pd
import urllib.parse
import json
from escapejson import escapejson
import pyppeteer

HARVARD_LOGIN_URL = "https://www.pin1.harvard.edu/cas/login"
ORDERS_URL = "https://sysbiolabops.hms.harvard.edu/orders/paulsson"
AUTH_BUTTON_SELECTOR = "div.push-label button.auth-button"
LAST_PAGE_SELECTOR = "li.pager__item--last > a"
ORDERS_TABLE_SELECTOR = "div#block-warplab1-content div.view-content table.views-table"


async def _get_table(browser, url, selector="table"):
    page = await browser.newPage()
    await page.goto(url)
    await page.waitForSelector(selector)
    table = await page.querySelectorEval(selector, "(element) => element.outerHTML")
    return pd.read_html(table)[0]


# SEE: https://beenje.github.io/blog/posts/parsing-javascript-rendered-pages-in-python-with-pyppeteer/
async def get_orders(
    browser,
    created_min=None,
    created_max=None,
    description=None,
    vendor=None,
    catalog=None,
    price_min=None,
    submitted=None,
    notes=None,
):
    page = await browser.newPage()
    params = {
        "created[min]": created_min,
        "created[max]": created_max,
        "field_description_value": description,
        "field_vendors_target_id": vendor,
        "field_total_price_value": price_min,
        "combine": submitted,
        "field_notes_value": notes,
    }
    for k, v in list(params.items()):
        if v is None:
            del params[k]
    await page.goto(ORDERS_URL + "?" + urllib.parse.urlencode(params))
    num_pages = None
    try:
        last_page_url = await page.querySelectorEval(
            LAST_PAGE_SELECTOR, "(element) => element.href", timeout=500
        )
    except:
        num_pages = 1
    if num_pages is None:
        num_pages = int(urllib.parse.parse_qs(last_page_url)["page"][0])
    param_sets = [params]
    # if there is one page, it breaks if we set page=1
    # so we leave the page parameter unset for the first page
    for page_num in range(2, num_pages + 1):
        page_params = params.copy()
        page_params["page"] = page_num
        param_sets.append(page_params)
    urls = [ORDERS_URL + "?" + urllib.parse.urlencode(params)]
    dfs = [await _get_table(browser, url, ORDERS_TABLE_SELECTOR) for url in urls]
    df = pd.concat(dfs, ignore_index=True)
    df = df.iloc[:, 1:]
    for col in ("Submitted", "Ordered"):
        df[col] = pd.to_datetime(df[col])
    return df


async def login_harvard(username, password, service, browser=None):
    if not browser:
        browser = await pyppeteer.launch()
    page = await browser.newPage()
    login_url = HARVARD_LOGIN_URL + f"?service={urllib.parse.quote(service)}"
    await page.goto(login_url)
    await page.waitForSelector("input#username")
    # escape
    username = escapejson(json.dumps(username))
    password = escapejson(json.dumps(password))
    login_script = f"""
    () =>
    {{
        document.getElementById('username').value = {username};
        document.getElementById('password').value = {password};
        document.getElementById('submitLogin').click();
    }}
    """
    await page.evaluate(login_script)
    await page.waitForSelector("iframe#duo_iframe")
    duo_iframe = next(frame for frame in page.frames if frame.name == "duo_iframe")
    await duo_iframe.waitForSelector(AUTH_BUTTON_SELECTOR, visible=True)
    auth_button = await duo_iframe.querySelector(AUTH_BUTTON_SELECTOR)
    await auth_button.click()
    return browser
