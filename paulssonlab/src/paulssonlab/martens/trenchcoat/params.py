from ruamel.yaml import YAML
import tables
import pathlib


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
