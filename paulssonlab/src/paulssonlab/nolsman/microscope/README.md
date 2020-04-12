# Microscope automation

Code for automating image acquisition with MicroManager's Python API.

## Installation

By default, `.envrc` sets MICROMANAGER_DIR to `/c/Program Files/Micro-Manager-2.0gamma`. If necessary, change MICROMANAGER_DIR to the appropriate value, ensuring that the path is in Cygwin style (beginning with `/c/` instead of `C:\`) and is surrounded by quotes (and no backslashes to escape spaces). After editing `.envrc`, you will have to run `direnv allow`.

Note that the `paulssonlab` git commit hooks do not work on Python 2, so you will have to `conda activate paulssonlab` before you `git commit` (this can also be accomplished by `cd $root`). This will be remedied soon, once MicroManager becomes Python 3-compatible.

## Contributors

- Noah Olsman
