import pandas as pd
import io
from requests_html import HTMLSession


def get_semrock_spectra(urls, bins=None):
    if isinstance(urls, str):
        urls = [urls]
    session = HTMLSession()
    spectra = {}
    for url in urls:
        html = session.get(url).html
        catalog_numbers = [
            a.text for a in html.find("#resultsView .cartSection > h1 > a")
        ]
        for catalog_number in catalog_numbers:
            spectrum_number = "-".join(catalog_number.split("-")[:2]).replace("/", "_")
            name = f"Semrock {spectrum_number}"
            spectrum_urls = [
                f"https://www.semrock.com/_ProductData/Spectra/{spectrum_number}_Spectrum.txt",
                f"https://www.semrock.com/_ProductData/Spectra/{spectrum_number}_DesignSpectrum.txt",
            ]
            spectrum = None
            for spectrum_url in spectrum_urls:
                res = session.get(spectrum_url)
                if not res.ok:
                    continue
                lines = io.StringIO(
                    "\n".join(
                        [l for l in res.text.split("\n") if not l.startswith("---")]
                    )
                )
                spectrum = pd.read_csv(
                    lines,
                    sep="\t",
                    skiprows=4,
                    names=["wavelength", spectrum_number],
                    index_col=0,
                ).squeeze()
                break
            if spectrum is None:
                raise ValueError(f"could not find spectrum for '{spectrum_number}'")
            if bins is not None:
                spectrum = interpolate_dataframe(spectrum, bins)
            spectra[name] = spectrum
    return spectra
