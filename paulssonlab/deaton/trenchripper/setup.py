import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="TrenchRipper",
    version="0.2.5",
    author="Daniel Eaton",
    packages=setuptools.find_packages(),
    install_requires=[
        "pandas",
        "numpy",
        "bokeh",
        "h5py",
        "h5py-cache==1.0",
        "scipy",
        "scikit-image",
        "jupyter-client",
        "jupyter-core",
        "jupyterlab",
        "jupyterlab-server",
        "jupyter-server-proxy==1.5.0",
        "matplotlib",
        "dask[complete]" "dask-jobqueue",
        "tifffile",
        "ipywidgets",
        "pulp",
        "fastparquet",
        "torch",
        "torchvision",
        "tables",
        "scikit-learn",
        "seaborn",
        "nd2reader==3.2.3",
        "parse",
        "qgrid",
        "holoviews[recommended]",
        "opencv-python",
    ],
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/DanielScottEaton/TrenchRipper",
    python_requires=">=3.7",
)

# BROKEN
#     - bokeh==2.0.2
#     - chardet==3.0.4
#     - click==7.1.2
#     - cloudpickle==1.4.0
#     - distributed==2.15.1
#     - idna==2.9
#     - pexpect==4.8.0
#     - pytz==2020.1
#     - six==1.14.0
#     - typing-extensions==3.7.4.2

# WORKING
#     - bokeh==1.4.0
#     - click==7.1.1
#     - cloudpickle==1.3.0
#     - distributed==2.12.0
#     - importlib-metadata==1.6.0
#     - pytz==2019.3
#     - zipp==3.1.0
