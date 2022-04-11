import pandas as pd
import holoviews as hv
import param
import requests
from requests_html import HTMLSession
from IPython.display import display
import io
import re
import paulssonlab.shenker.microscope_design.util as util


def get_fpbase_spectrum(id_, bins=None):
    df = pd.read_csv(
        f"https://www.fpbase.org/spectra_csv/?q={id_}", index_col="wavelength"
    ).squeeze()
    if bins is not None:
        df = util.interpolate_dataframe(df, bins)
    return df


def get_fpbase_protein_spectra(bins=None):
    fpbase_entries = requests.get("https://www.fpbase.org/api/proteins/spectra/").json()
    fps = {}
    for fpbase_entry in fpbase_entries:
        fp = {}
        spectra = []
        for spectrum_type in fpbase_entry["spectra"]:
            state_name = re.sub(r"^default_", "", spectrum_type["state"])
            for k, v in spectrum_type.items():
                if k == "state":
                    continue
                elif k == "data":
                    spectra.append(
                        pd.DataFrame(v, columns=["wavelength", state_name]).set_index(
                            "wavelength"
                        )
                    )
                else:
                    fp[f"{state_name}_{k}"] = v
        df = pd.concat(spectra, axis=1)
        if bins is not None:
            df = util.interpolate_dataframe(df, bins)
        fp["spectra"] = df
        fps[fpbase_entry["name"]] = fp
    return fps


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
                spectrum = util.interpolate_dataframe(spectrum, bins)
            spectra[name] = spectrum
    return spectra


class SpectraViewer(param.Parameterized):
    fp = param.Selector(label="fluorophore")
    dc = param.Selector(label="dichroic")
    lp = param.Selector(label="longpass")
    x_range = (None, None)
    y_range = (None, None)

    def __init__(
        self,
        fps,
        dichroics,
        longpass_filters,
        fp_names=None,
        dc_names=None,
        lp_names=None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.fps = fps
        self.dichroics = dichroics
        self.longpass_filters = longpass_filters
        if fp_names is None:
            fp_names = fps.keys()
        if dc_names is None:
            dc_names = dichroics.keys()
        if lp_names is None:
            lp_names = longpass_filters.keys()
        self.param["fp"].objects = list(fp_names)
        self.param["dc"].objects = list(dc_names)
        self.param["lp"].objects = list(lp_names)
        for k in ("fp", "dc", "lp"):
            setattr(self, k, self.param[k].objects[0])

    def keep_zoom(self, x_range, y_range):
        self.x_range = x_range
        self.y_range = y_range

    @param.depends("fp", "dc", "lp")
    def view(self):
        # plot = (
        #     self.fps[self.fp]["spectra"].fillna(0).hvplot()
        #     * self.dichroics[self.dc].hvplot()
        #     * self.longpass_filters[self.lp].hvplot()
        # )
        plot = (
            hv.Curve(self.fps[self.fp]["spectra"]["ex"].fillna(0).values).opts(
                color=hv.Cycle.default_cycles["default_colors"][0]
            )
            * hv.Curve(self.fps[self.fp]["spectra"]["em"].fillna(0).values).opts(
                color=hv.Cycle.default_cycles["default_colors"][1]
            )
            * hv.Curve(self.dichroics[self.dc].values).opts(
                color=hv.Cycle.default_cycles["default_colors"][2]
            )
            * hv.Curve(self.longpass_filters[self.lp].values).opts(
                color=hv.Cycle.default_cycles["default_colors"][3]
            )
        )
        plot = plot.opts(height=400, width=700)
        # plot = plot.opts(hv.opts.Curve(tools=['vline'])) # TODO: need to figure out a way to get all lines in the same tooltip... maybe make a dataframe and a single Curve object?
        plot = plot.redim.range(x=self.x_range, y=self.y_range)
        # FROM: https://discourse.holoviz.org/t/keep-zoom-level-when-changing-between-variables-in-a-scatter-plot/1120/2
        rangexy = hv.streams.RangeXY(
            source=plot, x_range=self.x_range, y_range=self.y_range
        )
        rangexy.add_subscriber(self.keep_zoom)
        return plot


def show_heatmap(df, highlight_negative=False, **kwargs):
    with pd.option_context("display.max_rows", None, "display.max_columns", None):
        df = df.style.format(precision=2).background_gradient(
            **{"cmap": "RdPu", "axis": None, **kwargs}
        )
        if highlight_negative:

            def style_negative(v, props=""):
                return props if v < 0 else None

            df = df.applymap(
                style_negative,
                props=f"background-color:{highlight_negative};color:white;",
            )
        display(df)
