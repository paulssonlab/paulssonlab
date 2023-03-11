import pandas as pd
import re
from cytoolz import dissoc
import requests


def get_fpbase_proteins(**params):
    fpbase_entries = requests.get(
        "https://www.fpbase.org/api/proteins/", params={**params, "format": "json"}
    ).json()
    fps = {}
    for fpbase_entry in fpbase_entries:
        fp = dissoc(fpbase_entry, "name")
        fp["states"] = {
            state["name"]: dissoc(state, "name") for state in fpbase_entry["states"]
        }
        fps[fpbase_entry["name"]] = fp
    return fps


def get_fpbase_spectrum(id_):
    df = pd.read_csv(
        f"https://www.fpbase.org/spectra_csv/?q={id_}", index_col="wavelength"
    ).squeeze()
    return df


def get_fpbase_protein_spectra():
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
        fp["spectra"] = df
        fps[fpbase_entry["name"]] = fp
    return fps
