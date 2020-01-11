# Paulsson Lab
This repo contains code and notebooks used in the Paulsson Lab.

## Installation
After cloning this repository, run the following (replace `paulssonlab` with the desired name for your conda environment):
```
conda env create -n paulssonlab -f environment.yml
echo "conda activate paulssonlab" > .envrc
direnv allow
pre-commit install
nbstripout --install
```

The first command creates a conda environment using the list of packages in `environment.yml`. The next two commands make this conda environment activate automatically when you `cd` into this repo (if you aren't using direnv, skip these two commands and run an explicit `conda activate paulssonlab` instead). The last two commands set up the pre-commit hooks that automatically format code (in both `.ipynb` and `.py` files) and strip output from jupyter notebooks before adding them to git. If you're starting from an empty conda environment, install these tools with `conda install -c conda-forge pre-commit black jupytext nbstripout`.
