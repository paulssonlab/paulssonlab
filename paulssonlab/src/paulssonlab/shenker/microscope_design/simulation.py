import numpy as np
import pandas as pd
import xarray as xr
import scipy
import pint

# FROM: https://pint.readthedocs.io/en/latest/tutorial.html#using-pint-in-your-projects
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity


def image_to_xarray(img, scale):
    xs = scale * np.arange(img.shape[1])
    ys = scale * np.arange(img.shape[0])[::-1]
    return xr.DataArray(img, coords=dict(x=xs, y=ys), dims=["y", "x"])


def offset_xarray(a, b, offsets):
    offsets = {name: getattr(b, name) + val for name, val in offsets.items()}
    return a.interp_like(b.assign_coords(**offsets)).assign_coords(b.coords)


def draw_excitation_line(
    width,
    edge_defocus,
    base_defocus,
    falloff,
    width_px=6500,
    height_px=300,
    height_sigma=3,
):
    if not ((0 <= falloff) and (falloff <= 1)):
        raise ValueError("falloff must be between 0 and 1")
    # expect defocus parameters in um
    width = float(width / ureg.um)
    edge_defocus = float(edge_defocus / ureg.um)
    base_defocus = float(base_defocus / ureg.um)
    x_dependence = np.abs(np.linspace(-1, 1, width_px)) ** 2
    sigma = edge_defocus * x_dependence + base_defocus
    x_max = width / 2
    xs = np.linspace(-x_max, x_max, width_px)
    y_max = height_sigma * sigma.max()
    ys = np.linspace(-y_max, y_max, height_px)
    #     img = scipy.stats.norm.pdf(
    #         np.arange(height_px)[:, np.newaxis], height_px / 2, sigma,
    #     ) * (1 - falloff * x_dependence)
    img = scipy.stats.norm.pdf(ys[:, np.newaxis], 0, sigma) * (
        1 - falloff * x_dependence
    )
    return xr.DataArray(img, coords=dict(x=xs, y=ys), dims=["y", "x"])


def bin_spectrum(df, bins):
    bin_assignment = pd.cut(df.index, bins).rename_categories(
        (bins.right + bins.left) / 2
    )
    return df.groupby(bin_assignment).mean()


def load_fpbase_spectra(urls, bins=None):
    spectra = {
        name: pd.read_csv(url)
        .rename(columns={f"{name} {kind}": kind for kind in ("ex", "em", "2p")})
        .set_index("wavelength")
        for name, url in urls.items()
    }
    if bins is not None:
        binned_spectra = {
            name: bin_spectrum(spectrum, bins) for name, spectrum in spectra.items()
        }
    return spectra


def clean_multiindex_csv(df):
    columns = pd.DataFrame(df.columns.tolist())
    columns.loc[columns[0].str.startswith("Unnamed:"), 0] = np.nan
    columns[0] = columns[0].fillna(method="ffill")
    mask = pd.isnull(columns[0])
    columns[0] = columns[0].fillna("")
    columns.loc[mask, [0, 1]] = columns.loc[mask, [1, 0]].values
    df.columns = pd.MultiIndex.from_tuples(columns.to_records(index=False).tolist())
    return df


def load_filter_spectra(filename):
    filters_all = clean_multiindex_csv(pd.read_csv(filename, header=[2, 3]).iloc[:, 1:])
    filters = {}
    filter_peaks = {}
    for pos in filters_all.columns.levels[0]:
        name = float(re.sub(r"(?:Filter|Sample 2E) @ ([\d.]+) mm", r"\1", pos))
        filters[name] = filters_all[pos].set_index(filters_all[pos].columns[0])
        filters[name].columns = ["transmission"]
        filters[name].index.name = "wavelength"
        filters[name] /= 100  # convert to fraction
        filter_peaks[name] = filters[name].iloc[:, 0].idxmax()
    return filters, filter_peaks
