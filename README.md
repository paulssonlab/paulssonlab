# Paulsson Lab
This repo contains code and notebooks used in the Paulsson Lab.

## Installation
You should fork this repo (click the button in the upper-right of https://github.com/paulssonlab/paulssonlab) under your personal github account. Then clone your personal fork with, e.g., `git clone git@github.com:shenker/paulssonlab.git`. After cloning this repository, run the following (replace `paulssonlab` with the desired name for your conda environment):
```
cd paulssonlab
git remote add upstream git@github.com:paulssonlab/paulssonlab.git
conda env create -n paulssonlab -f environment.yml
echo "conda activate paulssonlab" > .envrc
direnv allow
pre-commit install
pre-commit install -t commit-msg
pre-commit install-hooks
cd ..
git clone git@github.com:shenker/nbcleanse.git
cd paulssonlab
python ../nbcleanse/nbcleanse.py install
```

First we `cd` into the repo directory. The second command allows you to push to the central paulssonlab fork using `git push upstream` (use `git push origin` to push to your personal fork). The third command creates a conda environment using the list of packages in `environment.yml`. The next two commands make this conda environment activate automatically when you `cd` into this repo (if you aren't using direnv, skip these two commands and run an explicit `conda activate paulssonlab` instead). The remaining commands set up the pre-commit hooks that automatically format code (in both `.ipynb` and `.py` files) and strip output from jupyter notebooks before adding them to git.

## Basic Git workflow
When you make changes, `git add path/to/modified/file` to add them to the git index. When you have added all related changes, `git commit` (or `g ci`) them. Follow [these best practices](https://chris.beams.io/posts/git-commit/) for writing informative git commit messages. To push to your own fork, `git push origin` (by default, `origin` is the default remote, so you can just `git push`).

When you want to pull the latest changes from the rest of the lab, `git pull upstream master`. When you are ready to share your changes with the rest of the lab, first `git pull upstream master` (and fix merge conflicts if any arise), then `git push upstream`.

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
