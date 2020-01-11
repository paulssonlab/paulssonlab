# template-repo
Use this as a starting point for any Paulsson Lab repo. Comes with code formatting pre-commit hooks (black, jupytext, nbstripout) preconfigured.

## Installation
After cloning this repository, run the following (replace `template` with the desired name for your conda environment):
```
conda env create -n template -f environment.yml
echo "conda activate template" > .envrc
direnv allow
pre-commit install
nbstripout --install
```

The first command creates a conda environment using the list of packages in `environment.yml`. The next two commands make this conda environment activate automatically when you `cd` into this repo (if you aren't using direnv, skip these two commands). The last two commands set up the pre-commit hooks that automatically format code (in both `.ipynb` and `.py` files) and strip output from jupyter notebooks before adding them to git. If you're starting from an empty conda environment, install these tools with `conda install -c conda-forge pre-commit black jupytext nbstripout`.
