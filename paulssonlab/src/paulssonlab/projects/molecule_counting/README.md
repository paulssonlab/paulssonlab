# molecule-counting
Counting molecules with epifluorescence microscopy using photobleaching and binomial partitioning.

## Installation
```
initenv molecule_counting
jupyter labextension install @jupyter-widgets/jupyterlab-manager @pyviz/jupyterlab_pyviz
direnv allow
echo `git rev-parse --show-toplevel`/paulssonlab/src > `python -c 'import site; print(site.getsitepackages()[0])'`/paulssonlab.pth
```

## Contributors

- Noah Olsman
- Jacob Quinn Shenker
