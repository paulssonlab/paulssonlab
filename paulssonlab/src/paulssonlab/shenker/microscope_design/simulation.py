import numpy as np
import pandas as pd


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
