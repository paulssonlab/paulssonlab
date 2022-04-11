# !/usr/bin/env python

from setuptools import find_packages, setup

setup(
    name='paulssonlab',
    version='0.0.1',
    description='Analysis tools for the Paulsson Lab at Harvard Medical School',
    author='Jacob Shenker',
    author_email='aroy@alum.mit.edu',
    url='https://github.com/paulssonlab/paulssonlab',
    python_requires='>=3.8',
    packages=find_packages(),
    include_package_data=True,
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering",
    ],
)