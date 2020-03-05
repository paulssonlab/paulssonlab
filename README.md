# Paulsson Lab
This repo contains code and notebooks used in the Paulsson Lab.

## Installation
You should fork this repo (click the button in the upper-right of https://github.com/paulssonlab/paulssonlab) under your personal github account. Then clone your personal fork with, e.g., `git clone git@github.com:shenker/paulssonlab.git`. After cloning this repository, run the following (replace `paulssonlab` with the desired name for your conda environment):
```
cd paulssonlab
git remote add upstream git@github.com:paulssonlab/paulssonlab.git
pre-commit install
pre-commit install -t commit-msg
pre-commit install-hooks
cd ..
git clone git@github.com:shenker/nbcleanse.git
cd paulssonlab
conda env create -n nbcleanse -f nbcleanse-environment.yml
conda activate paulssonlab
python ../nbcleanse/nbcleanse.py install
conda deactivate
```

The first command allows you to push to the central paulssonlab fork using `git push upstream` (use `git push origin` to push to your personal fork). The remaining commands set up the pre-commit hooks that automatically format code (in both `.ipynb` and `.py` files) and strip output from jupyter notebooks before adding them to git.

## How to make a new project
TODO: How to structure python modules (example dir) so they can be easily imported. Primer on Python package structure.

## How to import an existing git repo
To import an existing git repo into the main `paulssonlab` monorepo (preserving commit history), first we rewrite the commit history to clean up Python and Jupyter files. Then we use `git-filter-repo` to rewrite history to move all files to a subdirectory. Then we merge this repo's commit history with this repo.
1. `conda activate nbcleanse` and install git-filter-repo with `pip install git-filter-repo`.
2. `cd path/to/nbcleanse` (where you cloned `nbcleanse` above)
3. `cd ..`
4. `git clone git@github.com:shenker/old-repo.git`
5. `cd old-repo`
6. Filter old-repo with `python ../nbcleanse/nbcleanse.py filter_repo` (this will take a few minutes).
7. Run `git filter-repo --strip-blobs-bigger-than 2M --to-subdirectory-filter shenker/old-repo`
8. Then merge this repo:
```
cd path/to/paulssonlab # this repo
git remote add -f old-repo path/to/old-repo
git merge --no-verify --allow-unrelated-histories old-repo/master
git remote rm old-repo
```
