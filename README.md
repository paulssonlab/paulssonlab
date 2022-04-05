# Paulsson Lab
This repo contains code and notebooks used in the Paulsson Lab.

## Installation
Open a shell, `cd` to the directory where you want to clone the repo, and run the following
```
/bin/bash -c "$(curl -fsSL https://gist.githubusercontent.com/shenker/11cccba44d843c5c135a152f55ef9d51/raw/paulssonlab-install.sh)"
```

If you are running this on O2, you *must* do this in an interactive job. You can start one with `irun1`. If you have not configured your O2 account with the Paulsson lab [dotfiles](https://github.com/paulssonlab/dotfiles), do that first.

Note that you need a `bash`-compatible shell installed (this is a given on Linux/macOS systems, but on Windows you will need to use WSL or scoop to set up `bash`).

## How to make a new project
Science projects for which you are the only contributor should be kept in a submodule named after your last name (optionally with your first initial as well), e.g., `paulsson/src/paulssonlab/shenker/my_project` For shared science projects that involve multiple contributors, notebooks and code specific to that science project should kept in a submodule under the `paulssonlab.projects` hierarchy, e.g., `paulsson/src/paulssonlab/projects/my_project`. Infrastructure code (segmentation, tracking, etc.) that is widely useful and not specific to any one science project should be kept in the base `paulssonlab` hierarchy, e.g., `paulssonlab/src/paulssonlab/segmentation`.

To make a new project, cd to an empty directory in an appropriate location and run `make-project`. For example, to make a new personal science project, you could run:
```
mkdir $src/shenker/my_project
cd $src/shenker/my_project
make-project
```
Here we the convenient fact that `$src` is defined to be `path/to/repo/paulssonlab/src/paulssonlab` (see below).

 Each project should contain a `README.md` file describing the project (summary, contributors, literature references, and any special installation instructions), an `environment.yml` file containing the `conda` packages required by that project, and an `.envrc` file that causes the conda environment to be created/activated. `make-project` will create these three files using default values. In particular, you will need to edit the `README.md`. By default, the `.envrc` contains two lines:
```
source_up .envrc_base
initenv
```
The first line must be `source_up .envrc_base`, this ensures that:
- `$root` is defined to be the path to the root of the git repo
- `$src` is defined to be `$root/paulssonlab/src/paulssonlab`
- `initenv` is defined
The second line, calling `initenv` with no arguments, will ask the user for the name of the conda environment to create. If an environment with that name already exists, it will be activated; if not, the user will be asked if they want to create it. The resulting environment will contain all the packages specified in `environment.yml`. The name of the environment will be written to `.envname`, so that the next time the user `cd`s to this directory, that environment will be activated automatically. By default, `initenv` uses `environment.yml` file in the same directory as the `.envrc` file, but you can specify an alternative file as the first argument to `initenv`: e.g., `initenv path/to/environment.yml`.

If you let `initenv` create all your conda environments, the top-level `paulssonlab` package will automatically be added to the `PYTHONPATH`. This allows you to easily import code from any submodule of the `paulssonlab` package, not just code within this project directory. It is strongly suggested that in your code (both Jupyter notebooks and `.py` files) you exclusively use absolute imports: e.g., `from paulssonlab.projects.my_project.segmentation import segment` instead of `from segmentation import segment`; many things will break in unexpected ways if you do not do this.

## Basic Git workflow
When you make changes, `git add path/to/modified/file` to add them to the git index. When you have added all related changes, `git commit` them. Follow [these best practices](https://chris.beams.io/posts/git-commit/) for writing informative git commit messages. When you commit, git runs pre-commit hooks to reformat code and check line endings. If one of the pre-commit hooks fails, it will modify the file on disk and abort the commit. If that happens, `git add path/to/modified/file` and try committing again. To push to your own personal fork, `git push origin`.

When you want to pull the latest changes from the rest of the lab, `git pull upstream`. When you are ready to share your changes with the rest of the lab, first `git pull upstream` (and fix merge conflicts if any arise), then `git push upstream`.

If you want to incorporate a change that a user (e.g., `nolsman`) has pushed to their own personal fork but not yet pushed to the main `paulssonlab` fork, you can `hub fetch nolsman`, then `git merge nolsman/master`.

If you push changes to your personal fork, you may want to pull them from another clone (e.g., if you have your fork cloned on both your laptop and O2). You can do this with `git pull origin master`.

## Common problems

- The `bioconda` channel (if you are using it) must be listed below `conda-forge` in `environment.yml` files, or else you will get errors about package conflicts.

## How to import an existing Git repo
To import an existing git repo into the main `paulssonlab` monorepo (preserving commit history), first we rewrite the commit history to clean up Python and Jupyter files. Then we use `git-filter-repo` to rewrite history to move all files to a subdirectory. Then we merge this repo's commit history with this repo.
1. `cd path/to/paulssonlab/.nbcleanse`
2. The `nbcleanse` environment has been created automatically for you; activate it with `conda activate nbcleanse`. If you need to create it, run `conda env create -n nbcleanse -f environment.yml` before activating it.
3. Ensure that the `nbcleanse` environment is activated and run `mamba install -c conda-forge git-filter-repo`
4. `cd ../..`
5. `git clone git@github.com:shenker/old-repo.git`
6. `cd old-repo`
7. Filter old-repo with `python ../paulssonlab/.nbcleanse/nbcleanse.py filter-repo` (this will take a few minutes).
8. Run `git filter-repo --strip-blobs-bigger-than 2M --to-subdirectory-filter shenker/old-repo`
9. Then merge this repo:
```
cd path/to/paulssonlab # this repo
git remote add -f old-repo path/to/old-repo
git merge --no-verify --allow-unrelated-histories old-repo/master
git remote rm old-repo
```
