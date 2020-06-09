from ruamel.yaml import YAML
import tables
import pathlib

# Input path to parameters YAML file, return dictionary of parameters
def read_params_file(file_path):
    yaml = YAML(typ="safe", pure=True)
    return yaml.load(pathlib.Path(file_path))


# Write dictionary of parameters to a YAML file
def write_params(params, file_path):
    yaml = YAML(typ="safe", pure=True)
    yaml.dump(params, pathlib.Path(file_path))
