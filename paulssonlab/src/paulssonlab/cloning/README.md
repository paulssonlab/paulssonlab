# Cloning automation

Code for planning cloning as well as maintaining registries of strains, plasmid maps, and primers.

## Installation
If you want to use the conda environment called `cloning`, cd into this directory and run the following.
```
initenv cloning
direnv allow
echo `git rev-parse --show-toplevel`/paulssonlab/src > `python -c 'import site; print(site.getsitepackages()[0])'`/paulssonlab.pth
```

You will also need to authenticate with the `paulssonlab@gmail.com` Google account to access Google Calendar, Google Sheets, and Google Drive. Ask Jacob how to do this.

## Contributors

- Jacob Quinn Shenker
