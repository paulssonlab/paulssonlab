import pathlib

import tables
from ruamel.yaml import YAML


def read_params_file(file_path):
    """
    Input path to parameters YAML file, return dictionary of parameters
    """
    yaml = YAML(typ="safe", pure=True)
    return yaml.load(pathlib.Path(file_path))


def write_params(params, file_path):
    """
    Write dictionary of parameters to a YAML file
    """
    yaml = YAML(typ="safe", pure=True)
    yaml.dump(params, pathlib.Path(file_path))


def read_params_string(string):
    """
    Input string with YAML data, return a dictionary of parameters
    """
    yaml = YAML(typ="safe", pure=True)
    return yaml.load(string)
