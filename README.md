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

## Project setup
Science projects that only you are working on should be kept in a submodule named after your last name (optionally with your first initial as well), e.g., `paulsson/src/paulssonlab/shenker/my_project` For shared science projects that more than one person are working on, notebooks and code specific to that science project should kept in a submodule under the `paulssonlab.projects` hierarchy, e.g., `paulsson/src/paulssonlab/projects/my_project`. Infrastructure code (segmentation, tracking, etc.) that is widely useful and not specific to any one science project should be kept in the base `paulssonlab` hierarchy, e.g., `paulssonlab/src/paulssonlab/segmentation`.

After setting up the repo (above), you will need to do the following setup for each science project you want to work on. `cd` to the project directory using, e.g., `cd paulsson/src/paulssonlab/projects/my_project`. Each project should contain a `README.md` file describing the project (summary, contributors, literature references, and any special installation instructions) and an `environment.yml` file containing the `conda` packages required by that project. From within this project directory, run the following:
```
conda env create -n my_project -f environment.yml
echo "conda activate my_project" > .envrc
# the following command only applies if you're using one of: ipywidgets, holoviews, or dask
jupyter labextension install @jupyter-widgets/jupyterlab-manager @pyviz/jupyterlab_pyviz
direnv allow
echo `git rev-parse --show-toplevel`/paulssonlab/src > `python -c 'import site; print(site.getsitepackages()[0])'`/paulssonlab.pth
```

The first line creates a new conda environment called `my_project` that contains all the conda packages listed in `environment.yml`. The next two lines ensure that this conda environment is loaded whenever we `cd` into this project directory. The last line adds the `paulssonlab` package to the PYTHONPATH within this conda environment. This allows you to easily import code from any submodule of the `paulssonlab` package, not just code within this project directory. It is strongly suggested that in your code (both Jupyter notebooks and `.py` files) you exclusively use absolute imports: e.g., `from paulssonlab.projects.my_project.segmentation import segment` instead of `from segmentation import segment`; many things will break in unexpected ways if you do not do this.

## Basic Git workflow
When you make changes, `git add path/to/modified/file` to add them to the git index. When you have added all related changes, `git commit` them. Follow [these best practices](https://chris.beams.io/posts/git-commit/) for writing informative git commit messages. To push to your own fork, `git push origin` (by default, `origin` is the default remote, so you can just `git push`).

When you want to pull the latest changes from the rest of the lab, `git pull upstream master`. When you are ready to share your changes with the rest of the lab, first `git pull upstream master` (and fix merge conflicts if any arise), then `git push upstream`. Alternatively, `hub sync` is a shortcut that does both `git pull origin` and `git pull upstream master`.

If you want to incorporate a change that a user (e.g., `nolsman`) has pushed to their own fork but not yet pushed to the main `paulssonlab` fork, you can `hub fetch nolsman`, then `git merge nolsman/master`.

## How to make a new project
TODO: How to structure python modules (example dir) so they can be easily imported. Primer on Python package structure.

## How to import an existing git repo
To import an existing git repo into the main `paulssonlab` monorepo (preserving commit history), first we rewrite the commit history to clean up Python and Jupyter files. Then we use `git-filter-repo` to rewrite history to move all files to a subdirectory. Then we merge this repo's commit history with this repo.
1. `cd path/to/nbcleanse` (where you cloned `nbcleanse` above)
2. `conda env create -n nbcleanse -f environment.yml` (or if you have already created the `nbcleanse` conda environment, you can `conda activate nbcleanse`)
3. `cd ..`
4. `git clone git@github.com:shenker/old-repo.git`
5. `cd old-repo`
6. Filter old-repo with `python ../nbcleanse/nbcleanse.py filter-repo` (this will take a few minutes).
7. Run `git filter-repo --strip-blobs-bigger-than 2M --to-subdirectory-filter shenker/old-repo`
8. Then merge this repo:
```
cd path/to/paulssonlab # this repo
git remote add -f old-repo path/to/old-repo
git merge --no-verify --allow-unrelated-histories old-repo/master
git remote rm old-repo
```
