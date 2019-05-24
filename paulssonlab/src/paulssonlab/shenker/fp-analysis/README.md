# fp-analysis
Data analysis scripts for the fluorescent protein library project.

## Installation
After cloning this repository, run the following:
```
conda env create -n fp-analysis -f environment.yml
echo "conda activate fp-analysis" > .envrc
direnv allow
pre-commit install
nbstripout --install
```