# Microscope automation

Code for automating image acquisition with MicroManager's Python API.

## Installation
If you want to use the conda environment called `deaton-microscope`, cd into this directory and run the following. Change MICROMANAGER_DIR to the appropriate value, ensuring that the path is in Cygwin style (beginning with `/c/` instead of `C:\`) and is surrounded by quotes (and no backslashes to escape spaces).
```
initenv deaton-microscope
echo 'export NOTEBOOK_DIR=$PWD' >> .envrc
echo 'export MICROMANAGER_DIR="/c/Program Files/Micro-Manager-2.0gamma"' >> .envrc
direnv allow
echo `git rev-parse --show-toplevel`/paulssonlab/src > `python -c 'import site; print(site.getsitepackages()[0])'`/paulssonlab.pth
```

Note that the `paulssonlab` git commit hooks do not work on Python 2, so you will have to `conda activate paulssonlab` before you `git commit`.

## Contributors

- Daniel Eaton
