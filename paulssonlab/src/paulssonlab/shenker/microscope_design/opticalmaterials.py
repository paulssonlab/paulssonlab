# import numpy as np
import jax
from jax import numpy as np
import pandas as pd

# from cytoolz import partial
from jax.tree_util import Partial
from scipy.interpolate import InterpolatedUnivariateSpline
import zipfile
import yaml
import os
import io

# NOTE: wavelengths are in µm!

# https://github.com/quartiq/rayopt/blob/master/rayopt/rii.py
# and
# https://github.com/quartiq/rayopt/blob/master/rayopt/material.py
# were used as reference for the refractiveindex.info database format
# and refractive index formulae

# TODO: needs testing
def n_schott(c, w):
    n = c[0] + c[1] * w**2
    for i, ci in enumerate(c[2:]):
        n += ci * w ** (-2 * (i + 1))
    return np.sqrt(n)


# TODO: needs testing
def n_sellmeier(c, w):
    ndim = np.ndim(w)
    w = np.expand_dims(w, axis=-1)
    w2 = w**2
    c0, c1 = c.reshape([-1] + [1] * ndim + [2]).T
    return np.sqrt(1.0 + (c0 * w2 / (w2 - c1**2)).sum(axis=ndim))


# TODO: needs testing
def n_sellmeier_squared(c, w):
    ndim = np.ndim(w)
    w = np.expand_dims(w, axis=-1)
    w2 = w**2
    c0, c1 = c.reshape([-1] + [1] * ndim + [2]).T
    return np.sqrt(1.0 + (c0 * w2 / (w2 - c1)).sum(axis=ndim))


# TODO: needs testing
def n_sellmeier_squared_transposed(c, w):
    ndim = np.ndim(w)
    w = np.expand_dims(w, axis=-1)
    w2 = w**2
    c0, c1 = c.reshape([-1] + [1] * ndim + [2])
    return np.sqrt(1.0 + (c0 * w2 / (w2 - c1)).sum(axis=ndim))


# TODO: needs testing
def n_conrady(c, w):
    return c[0] + c[1] / w + c[2] / w**3.5


# TODO: needs testing
def n_herzberger(c, w):
    l = 1.0 / (w**2 - 0.028)
    return (
        c[0] + c[1] * l + c[2] * l**2 + c[3] * w**2 + c[4] * w**4 + c[5] * w**6
    )


def n_sellmeier_offset(c, w):
    ndim = np.ndim(w)
    w = np.expand_dims(w, axis=-1)
    w2 = w**2
    c0, c1 = c[1 : 1 + (c.shape[0] - 1) // 2 * 2].reshape([-1] + [1] * ndim + [2]).T
    return np.sqrt(1.0 + c[0] + (c0 * w2 / (w2 - c1**2)).sum(axis=ndim))


def n_sellmeier_squared_offset(c, w):
    ndim = np.ndim(w)
    w = np.expand_dims(w, axis=-1)
    w2 = w**2
    c0, c1 = c[1 : 1 + (c.shape[0] - 1) // 2 * 2].reshape([-1] + [1] * ndim + [2]).T
    return np.sqrt(1.0 + c[0] + (c0 * w2 / (w2 - c1)).sum(axis=ndim))


# TODO: needs testing
def n_handbook_of_optics1(c, w):
    return np.sqrt(c[0] + (c[1] / (w**2 - c[2])) - (c[3] * w**2))


# TODO: needs testing
def n_handbook_of_optics2(c, w):
    return np.sqrt(c[0] + (c[1] * w**2 / (w**2 - c[2])) - (c[3] * w**2))


# TODO: needs testing
def n_extended2(c, w):
    n = c[0] + c[1] * w**2 + c[6] * w**4 + c[7] * w**6
    for i, ci in enumerate(c[2:6]):
        n += ci * w ** (-2 * (i + 1))
    return np.sqrt(n)


# TODO: needs testing
def n_hikari(c, w):
    n = c[0] + c[1] * w**2 + c[2] * w**4
    for i, ci in enumerate(c[3:]):
        n += ci * w ** (-2 * (i + 1))
    return np.sqrt(n)


# TODO: needs testing
def n_gas(c, w):
    ndim = np.ndim(w)
    w = np.expand_dims(w, axis=-1)
    c0, c1 = c.reshape([-1] + [1] * ndim + [2])
    return 1.0 + (c0 / (c1 - w**-2)).sum(axis=ndim)


# TODO: needs testing
def n_gas_offset(c, w):
    return c[0] + n_gas(c[1:], w)


# TODO: needs testing
def n_refractiveindex_info(c, w):
    ndim = np.ndim(w)
    w = np.expand_dims(w, axis=-1)
    c0, c1 = c[9:].reshape([-1] + [1] * ndim + [2]).T
    return np.sqrt(
        c[0]
        + c[1] * w ** c[2] / (w**2 - c[3] ** c[4])
        + c[5] * w ** c[6] / (w**2 - c[7] ** c[8])
        + (c0 * w**c1).sum(axis=ndim)
    )


# TODO: needs testing
def n_retro(c, w):
    w2 = w**2
    a = c[0] + c[1] * w2 / (w2 - c[2]) + c[3] * w2
    return np.sqrt(2 + 1 / (a - 1))


# TODO: needs testing
def n_cauchy(c, w):
    ndim = np.ndim(w)
    w = np.expand_dims(w, axis=-1)
    c0, c1 = c[1:].reshape([-1] + [1] * ndim + [2]).T
    return c[0] + (c0 * w**c1).sum(axis=ndim)


# TODO: needs testing
def n_polynomial(c, w):
    return np.sqrt(n_cauchy(c, w))


# TODO: needs testing
def n_exotic(c, w):
    return np.sqrt(
        c[0] + c[1] / (w**2 - c[2]) + c[3] * (w - c[4]) / ((w - c[4]) ** 2 + c[5])
    )


RII_FORMULAS = {
    "formula 1": n_sellmeier_offset,
    "formula 2": n_sellmeier_squared_offset,
    "formula 3": n_polynomial,
    "formula 4": n_refractiveindex_info,
    "formula 5": n_cauchy,
    "formula 6": n_gas_offset,
    "formula 7": n_herzberger,
    "formula 8": n_retro,
    "formula 9": n_exotic,
}


def _interpolate(x, y):
    func = InterpolatedUnivariateSpline(x, y)
    return func, func.derivative()


def get_n(material, nanometers=False):
    for formula_name, formula_func in RII_FORMULAS.items():
        data = material["DATA"].get(formula_name, None)
        if data is not None:
            n_microns_func = Partial(formula_func, data["coefficients"])
            if nanometers:
                n_func = lambda w: n_microns_func(w / 1e3)
            else:
                n_func = n_microns_func
            grad_n_func = np.vectorize(jax.grad(n_func))
            return n_func, grad_n_func
    data = material["DATA"].get("tabulated n", None)
    if data is None:
        data = material["DATA"].get("tabulated nk", None)
    if data is None:
        raise ValueError(
            f"material {material['name']} missing refractive index coefficients and tabulated data"
        )
    table = data["data"]
    wavelengths = table.index.values
    if nanometers:
        wavelengths *= 1e3
    return _interpolate(wavelengths, table["n"].values)


def get_k(material, nanometers=False):
    data = material["DATA"].get("tabulated k", None)
    if data is None:
        data = material["DATA"].get("tabulated nk", None)
    if data is None:
        raise ValueError(
            f"material {material['name']} missing extinction coefficient tabulated data"
        )
    table = data["data"]
    wavelengths = table.index.values
    if nanometers:
        wavelengths *= 1e3
    return _interpolate(wavelengths, table["k"].values)


def _parse_rii_material(material):
    all_data = {}
    for data in material["DATA"]:
        if data["type"].startswith("formula"):
            data = {
                k: np.fromstring(v, sep=" ") if k != "type" else v
                for k, v in data.items()
            }
        elif data["type"].startswith("tabulated"):
            columns = ["wavelength"]
            if data["type"] == "tabulated n":
                columns.append("n")
            elif data["type"] == "tabulated k":
                columns.append("k")
            elif data["type"] == "tabulated nk":
                columns.extend(["n", "k"])
            else:
                raise ValueError(f"unknown data type: {data['type']}")
            table = pd.read_csv(
                io.StringIO(data["data"]),
                sep=r"\s+",
                names=columns,
                index_col="wavelength",
            ).drop_duplicates()
            data = {**data, "data": table}
        all_data[data["type"]] = data
    return {**material, "DATA": all_data}


def parse_rii_catalog(filename):
    catalog = {}
    if os.path.isfile(filename):
        # assume filename is a zip file
        rii_zip = zipfile.ZipFile(filename, mode="r")
        open_ = rii_zip.open
        base_path = "database"
    else:
        open_ = open
        base_path = filename
    with open_(os.path.join(base_path, "library.yml"), "r") as f:
        rii_library = yaml.safe_load(f)
    for shelf in rii_library:
        for book in shelf["content"]:
            if "DIVIDER" in book:
                div = book["DIVIDER"]
                continue
            for page in book["content"]:
                if "DIVIDER" in page:
                    continue
                yml_filename = os.path.join(base_path, "data", page["data"])
                try:
                    with open_(yml_filename, "r") as f:
                        material = yaml.safe_load(f)
                except:
                    print(f"error, skipping {yml_filename}")
                    continue
                material["yml"] = page["data"]
                # stringifying probably unnecessary, except for PAGE
                # which is sometimes parsed as an integer
                material["BOOK"] = str(book["BOOK"])
                material["PAGE"] = str(page["PAGE"])
                material["name"] = str(page["name"])
                material["div"] = div
                material = _parse_rii_material(material)
                catalog[(material["BOOK"], material["PAGE"])] = material
    return catalog


def transmittance(wavelength, k, thickness, n=None, reflections=None, nanometers=False):
    if nanometers:
        raise NotImplementedError
    if reflections not in (None, "single", "multiple"):
        raise ValueError("reflections must be one of: None, single, multiple")
    if reflections == "multiple" and n is None:
        raise ValueError("n must be specified for multiple reflections")
    alpha = 4 * np.pi * k / wavelength
    tau = np.exp(-thickness * alpha)
    if reflections is None:
        T = tau
    else:
        n2 = complex(np.abs(n), k)
        R = np.abs((1 - n2) / (1 + n2))
        if reflections == "single":
            T = tau * (1 - R) ** 2
        else:
            T = tau * (1 - R) ** 2 / (1 - (R * tau) ** 2)
    return T
