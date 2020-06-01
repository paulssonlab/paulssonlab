import numpy as np
import pandas as pd
import xarray as xr
import scipy
import pint
import re

# FROM: https://pint.readthedocs.io/en/latest/tutorial.html#using-pint-in-your-projects
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

THORLABS_COLUMN_RE = r"(?:% )?Reflectance(?: \((?:((?:(?![()%,]).)+), )?((?:(?![()%,]).)+)\)|, (\S+) \(%\))"


def valid_range(ary):
    is_valid = ~np.isnan(ary)
    idx1 = is_valid.argmax()
    idx2 = len(ary) - is_valid[::-1].argmax()
    return idx1, idx2


def import_fpbase_spectrum(url):
    df = pd.read_csv(url, index_col="wavelength")
    name = re.sub(r"\s*ex|em|2p\s*", "", df.columns[0])
    df.columns = [c.replace(f"{name} ", "") for c in df.columns]
    return df, name


def import_fpbase_spectra(urls):
    spectra = {}
    for url in urls:
        spectrum, name = import_fpbase_spectrum(url)
        spectra[name] = spectrum
    return spectra


def bin_size(x):
    return (x[-1] - x[0]) / len(x)


def interpolate_dataframe(df, bins, union=False):
    old_index = df.index.astype(np.float)
    new_index = pd.Index(bins, name=df.index.name)
    union_index = new_index.union(old_index)
    new_df = df.reindex(index=union_index)
    new_df.interpolate(method="nearest", inplace=True)
    if not union:
        new_df = new_df.loc[new_index]
    return new_df


def bin_dataframe(df, bins):
    bin_assignment = pd.cut(df.index, bins).rename_categories(
        (bins.right + bins.left) / 2
    )
    return df.groupby(bin_assignment).mean()


def seesaw_spectrum(spectrum, amount):
    x = np.zeros(len(spectrum))
    idx1, idx2 = valid_range(spectrum)
    x[idx1:idx2] = np.linspace(-1, 1, idx2 - idx1)
    x = x.reshape((-1,) + (1,) * (np.ndim(spectrum) - 1))
    shape = -(x - 1) * (x + 1) * x * 3 * np.sqrt(3) / 4
    new_spectrum = (amount * shape + 1 / 2) * spectrum
    new_spectrum /= np.nanmax(new_spectrum, axis=0)
    return new_spectrum


def xarray_like(ary, data):
    return xr.DataArray(data, coords=ary.coords, dims=ary.dims)


def image_to_xarray(img, scale):
    xs = scale * np.arange(img.shape[1])
    ys = scale * np.arange(img.shape[0])[::-1]
    return xr.DataArray(img, coords=dict(x=xs, y=ys), dims=["y", "x"])


def shift_xarray(ary, shifts):
    new_coords = {name: getattr(ary, name) + val for name, val in shifts.items()}
    return ary.assign_coords(**new_coords)


def shift_and_interp(a, b, shifts, method="linear"):
    new_coords = {name: getattr(b, name) + val for name, val in shifts.items()}
    return a.interp_like(b.assign_coords(**new_coords), method=method).assign_coords(
        b.coords
    )


def generalized_normal_pdf(x, loc=0, scale=1, p=2, sigma=None):
    if loc != 0:
        x = x - loc
    if sigma is not None:
        scale = sigma / np.sqrt(scipy.special.gamma(3 / p) / scipy.special.gamma(1 / p))
    return scipy.stats.gennorm.pdf(x / scale, p) / scale


def draw_excitation_line(
    width,
    height,
    edge_defocus,
    falloff=0,
    p_vertical=2,
    p_horizontal=2,
    width_px=6500,
    height_px=300,
    height_padding_factor=3,
    dtype=np.float32,
):
    # expect defocus parameters in um
    width = (
        np.atleast_1d((width / ureg.um).to("dimensionless").magnitude)
        .reshape((-1, 1, 1))
        .astype(dtype)
    )
    height = (
        np.atleast_1d((height / ureg.um).to("dimensionless").magnitude)
        .reshape((-1, 1, 1))
        .astype(dtype)
    )
    edge_defocus = (
        np.atleast_1d((edge_defocus / ureg.um).to("dimensionless").magnitude)
        .reshape((-1, 1, 1))
        .astype(dtype)
    )
    falloff = np.atleast_1d(falloff).reshape((-1, 1, 1)).astype(dtype)
    p_vertical = np.atleast_1d(p_vertical).reshape((-1, 1, 1)).astype(dtype)
    p_horizontal = np.atleast_1d(p_horizontal).reshape((-1, 1, 1)).astype(dtype)
    if not (np.all(0 <= falloff) and np.all(falloff <= 1)):
        raise ValueError("falloff must be between 0 and 1")
    width_max = width.max()
    xs_normalized = np.linspace(-1, 1, width_px, dtype=dtype)
    xs = width_max / 2 * xs_normalized
    xs_scaled = width_max / width * xs_normalized[np.newaxis, np.newaxis, :]
    scale = edge_defocus * np.abs(xs_normalized) ** 2 + height / 2
    y_max = height_padding_factor * scale.max()
    ys = np.linspace(-y_max, y_max, height_px, dtype=dtype)
    falloff_profile = (1 - falloff) + falloff * generalized_normal_pdf(
        xs_scaled, p=p_horizontal
    ).astype(dtype)
    img = (
        generalized_normal_pdf(
            ys[np.newaxis, :, np.newaxis], scale=scale, p=p_vertical
        ).astype(dtype)
        * falloff_profile
    )
    if img.shape[0] == 1:
        return xr.DataArray(img[0], coords=dict(x=xs, y=ys), dims=["y", "x"])
    else:
        return xr.DataArray(img, coords=dict(x=xs, y=ys), dims=["ex", "y", "x"])


def clean_multiindex_csv(df):
    columns = pd.DataFrame(df.columns.tolist())
    columns.loc[columns[0].str.startswith("Unnamed:"), 0] = np.nan
    columns[0] = columns[0].fillna(method="ffill")
    mask = pd.isnull(columns[0])
    columns[0] = columns[0].fillna("")
    columns.loc[mask, [0, 1]] = columns.loc[mask, [1, 0]].values
    df.columns = pd.MultiIndex.from_tuples(columns.to_records(index=False).tolist())
    return df


def read_filter_spectra(filename):
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


def read_thorlabs(filename):
    first_row = pd.read_excel(filename, sheet_name=0, header=None, nrows=1)
    if (~first_row.isnull()).values.sum() == 1:
        skiprows = 1
        if first_row[first_row.T.first_valid_index()].values[0].startswith("Variation"):
            header = (0, 1)
        else:
            header = 0
    else:
        skiprows = 0
        header = 0
    print("skip", skiprows)
    sheets = pd.read_excel(filename, sheet_name=None, skiprows=skiprows, header=header)
    # equivalent to usecols=lambda x: "Unnamed" not in x,
    # but works for MultiIndex
    dfs = {}
    for sheet_name, sheet in sheets.items():
        # print(">>", sheet.columns)
        # print("$$", sheet)
        # print(sheet.columns.get_level_values(0).str.startswith("Unnamed"))
        print("$", sheet.columns.get_level_values(0))
        sheet = sheet.loc[
            :, ~sheet.columns.get_level_values(0).str.startswith("Unnamed")
        ]
        for col in sheet.columns:
            print(">", col)
            if (isinstance(col, tuple) and col[0].startswith("Wavelength (µm)")) or (
                not isinstance(col, tuple) and col.startswith("Wavelength (µm)")
            ):
                print("FIX")
                sheet[col].values[:] *= 1000
        is_index = sheet.columns.get_level_values(0).str.startswith("Wavelength")
        print("index", is_index)
        splits = np.concatenate((np.where(is_index)[0], (len(sheet.columns),)))
        print("s", splits)
        if len(splits) > 2:
            column_sets = np.vstack((splits[:-1], splits[1:]))
        else:
            column_sets = [splits]
        print(column_sets, column_sets[0])
        for column_set in column_sets:
            df = sheet.iloc[:, slice(*column_set)]
            # print(df)
            print("COLS", df.columns.values)
            if isinstance(df.columns, pd.MultiIndex):
                col_names = df.columns.get_level_values(1)
                col_names0 = df.columns.get_level_values(0)
                col_names = [col_names0[0], *col_names[1:]]
            else:
                col_names = df.columns.get_level_values(0)
            print("COL NAMES", col_names)
            assert col_names[0].startswith("Wavelength")
            col_names = ["Wavelength", *col_names[1:]]
            m = re.match(THORLABS_COLUMN_RE, col_names[1])
            if m and m.group(1):
                sheet_name = m.group(1)
            # split df
            print(df.columns)
            col_names = [
                re.sub(
                    r"^(?:% Reflectance|Reflectance \(%\))$",
                    "Reflectance",
                    re.sub(THORLABS_COLUMN_RE, "\\2\\3", c),
                )
                for c in col_names
            ]
            df.columns = pd.Index(col_names)
            # df.columns.set_levels(col_names, level=0, inplace=True)
            ###TODO
            print(df.columns.values)
            df.set_index("Wavelength", inplace=True)
            df = df.loc[: df.last_valid_index()]
            df.sort_index(inplace=True)
            df.values[:] /= 100  # % to fraction
            dfs[sheet_name] = df
    if len(dfs) == 1:
        merged_df = next(iter(dfs.values()))
    else:
        merged_df = pd.concat(dfs, axis=1)
    return merged_df
