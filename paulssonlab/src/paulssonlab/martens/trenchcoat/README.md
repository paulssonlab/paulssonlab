# TrenchCoat

## What is TrenchCoat?

TrenchCoat is a collection of command-line utilities which together process, analyze & display microscopy images. TrenchCoat is written in [Python](https://www.python.org/), uses [HDF5](http://www.pytables.org) for storage, [Napari](https://napari.org) for visualization, and supports multi-processing for improved performance.

## How do I install TrenchCoat?

Although TrenchCoat has mostly been tested on Linux and Windows, it should be cross-platform.

TrenchCoat requires the following Python packages:

  - click -- for generating the command-line interface
  - dask -- for displaying multi-dimensional datasets in Napari
  - ipywidgets -- for (optionally) enabling widget sliders on segmentation algorithms, for quick testing of parameters (may be removed in future versions)
  - magicgui -- for creating GUIs in Napari (work in progress)
  - napari -- interface for displaying microscopy images
  - pyside2 -- Qt library for Napari
  - numba -- used in the registration code (work in progress)
  - numpy -- for working with arrays
  - pathlib -- for handling paths
  - pytables -- for working with HDF5 files & data tables
  - ruamel.yaml -- for conveniently storing & loading text-based parameter files
  - scikit-image -- for analyzing images
  - scipy -- for analyzing images
  - tqdm -- for progress bars
  - xmltodict -- for converting some of the ND2 metadata
  - matplotlib -- for access to the tab20 colormap to display cell lineages

In addition, it requires a slightly modified version of the [nd2reader](https://github.com/jimrybarski/nd2reader) Python package (see [here](https://github.com/paulssonlab/paulssonlab/tree/master/paulssonlab/src/paulssonlab/martens/nd2reader)). See [below](#nd2reader) for an explanation.

The source can be downloaded from GitHub and executed from the command line like any Python program.

## How do I run TrenchCoat?

Run ``trenchcoat --help`` to display the following:

```
Usage: trenchcoat [OPTIONS] COMMAND [ARGS]...

  Invoke a Command to perform the desired operation:

Options:
  --help  Show this message and exit.

Commands:
  browse-hdf5               Use Napari to browse a dataset & to visualize...
  browse-kymographs         Use Napari to browse kymographs.
  browse-nd2                Use Napari to browse a directory of ND2 files.
  convert                   Convert a directory of ND2 files to an HDF5...
  correct-trench-numbering  Take into account when a trench lies near the...
  corrections               Generate camera bias and flat field corrections...
  kymographs                Generate kymographs.
  lineage-tracking          Track lineages using HDF5 tables.
  print-metadata            Print the HDF5 metadata to the terminal.
  segment                   Detect cells and write their properties and...
  trench-detect             Detect trenches and write their rectangular...
  trench-measurements       Analyze whole trenches, without cell...
```

Below is a description of each command, followed by an example series of steps used for a typical workflow. The commands have been ordered as one might perform them during analysis.

### print-metadata

```
Usage: trenchcoat print-metadata [OPTIONS]

  Print the HDF5 metadata to the terminal. If no sub-file is specified, then
  print a list of all sub-files.

Options:
  -i, --in-file TEXT        Input HDF5 file with metadata.  [default:
                            HDF5/data.h5; required]

  -f, --sub-file-name TEXT  Name of sub-file to print metadata for. This file
                            was a separate ND2 file before conversion to HDF5.

  --help                    Show this message and exit.
```

### browse-nd2

```
Usage: trenchcoat browse-nd2 [OPTIONS]

  Use Napari to browse a directory of ND2 files.

Options:
  -i, --in-dir TEXT               Input directory of ND2 files.  [default:
                                  ND2; required]

  -S, --napari-settings-file TEXT
                                  TrenchCoat settings file (YAML).  [default:
                                  trenchcoat_settings.yaml; required]

  --help                          Show this message and exit.
```

### convert

```
Usage: trenchcoat convert [OPTIONS]

  Convert a directory of ND2 files to an HDF5 file.

Options:
  -o, --out-dir TEXT           Output directory.  [default: HDF5; required]
  -i, --in-dir TEXT            Input directory of ND2 files.  [default: ND2;
                               required]

  -n, --num-cpu INTEGER RANGE  Number of CPUs to use.  [default: all CPUs]
  -f, --frames TEXT            List of frames to copy.  [default: all frames]
  -F, --fovs TEXT              List of FOVs to copy.  [default: all FOVs]
  --help                       Show this message and exit.
```

By default the entire ND2 datasets will be copied. However, conversion also supports only copying a subset of the ND2 datasets. This can be useful:
 - for testing purposes (reducing the time for analysis)
 - if it's known that only some of the data are needed

 For example, suppose one wanted to only include frames 1 and 5-7 (including 6), and FOVS 10-90 (inclusive):

```
trenchcoat convert -i ND2 -o HDF5 -f 1,5-7 -F 10-90
```

I most frequently use this feature for splitting up mother machine data into individual lanes, each withs its own FOV range. This helps to break up the data conceptually, and also comes in handy if the trench widths are different for the different lanes.

### browse-hdf5

```
Usage: trenchcoat browse-hdf5 [OPTIONS]

  Use Napari to browse a dataset & to visualize trenches and cell masks.

  camera_biases_file, flatfield_corrections_file are paths to HDF5 files
  containing channel-specific correction values.

  Add option to "compute" an image based on properties within segmented
  regions.

Options:
  -i, --images-file TEXT          Input HDF5 file with images.  [default:
                                  HDF5/; required]

  -m, --masks-file TEXT           Input HDF5 file with masks.
  -r, --regions-file TEXT         Input HDF5 file with regions.
  -S, --napari-settings-file TEXT
                                  TrenchCoat settings file (YAML).  [default:
                                  trenchcoat_settings.yaml; required]

  -V, --viewer-params TEXT        YAML file with image viewer params &
                                  information for displaying column data.

  -d, --data-table-file TEXT      Path to HDF5 file containing cell
                                  measurements for computed image.

  --help                          Show this message and exit.
```

### corrections

```
Usage: trenchcoat corrections [OPTIONS]

  Generate camera bias and flat field corrections matrices from images (work in progress).

Options:
  -o, --out-file TEXT             Output HDF5 file with corrections matrices.
                                  [default: corrections.h5; required]

  -i, --in-file TEXT              Input HDF5 file with images.  [default:
                                  BLANKS/data.h5; required]

  -D, --dark-channel TEXT         Name of dark channel (cannot be used with
                                  -B).

  -B, --background-values-file TEXT
                                  YAML file with background values (cannot be
                                  used with -D).

  --help                          Show this message and exit.
```

### trench-detect

```
Usage: trenchcoat trench-detect [OPTIONS]

  Detect trenches and write their rectangular regions to an HDF5 file.

Options:
  -o, --out-dir TEXT           Output directory.  [default: regions_table.h5; required]
  -i, --in-file TEXT           Input HDF5 file.  [default: HDF5/;
                               required]

  -n, --num-cpu INTEGER RANGE  Number of CPUs to use. [default: all CPUs]
  -P, --params-file TEXT       Regions detection parameters file (YAML).
                               [default: trench_params.yaml; required]

  --help                       Show this message and exit.
```

### segment

```
Usage: trenchcoat segment [OPTIONS]

  Detect cells and write their properties and masks to an HDF5 file.

Options:
  -o, --out-dir TEXT           Output for new HDF5 directory for masks &
                               measurements.  [default: SEGMENTATION;
                               required]

  -i, --in-file TEXT           Input HDF5 file with images.  [default:
                               HDF5/; required]

  -n, --num-cpu INTEGER RANGE  Number of CPUs to use. [default: all CPUs]
  -P, --params-file TEXT       Segmentation parameters file (YAML).  [default:
                               seg_params.yaml; required]

  -R, --regions-file TEXT      HDF5 file containing image regions. [default:
                               analyze entire image, no regions]

  --help                       Show this message and exit.
```

Note that *regions* are optional. Regions is a generic term for *e.g.* trenches. If no regions are specified, then the entire image is analyzed. This is useful for agar pads (which do not have trenches), or possibly for clustering-type trench detection which finds trenches after finding cells (currently unimplemented here).

### correct-trench-numbering

```
Usage: trenchcoat correct-trench-numbering [OPTIONS]

  Take into account when a trench lies near the edge of an image by re-
  numbering all the trenches as needed. This can be thought of as stage "drift"
  correction.

Options:
  -i, --in-file TEXT              Input HDF5 file with metadata.  [default:
                                  HDF5/; required]

  -R, --regions-file TEXT         HDF5 file containing image regions
                                  [default: regions_tables.h5]

  -O, --regions-out-file TEXT     Write renumbered trenches & regions.
                                  [required]

  -M, --measurements-file TEXT    Directory containing subdir with cell
                                  measurements with original trench numbering.
                                  [default: SEG/]

  -m, --measurements-out-file TEXT
                                  Write renumbered trenches paired with cell
                                  measurements.

  -H, --method TEXT               Method to use for calculating drift. Options
                                  are stage or image.  [default: stage;
                                  required]

  -C, --channel TEXT              Channel to use if using images to correct
                                  for drift.

  -T, --threshold INTEGER         Threshold for calling a renumbering event
                                  (pixels).  [default: 10; required]

  --help                          Show this message and exit.
```

### lineage-tracking

```
Usage: trenchcoat lineage-tracking [OPTIONS]

  Track lineages using HDF5 tables. Write a new HDF5 table.

Options:
  -i, --in-file TEXT           Input HDF5 file with cell measurements data and
                               trench labels.  [default: newsegs.h5; required]

  -o, --out-file TEXT          Output HDF5 file with lineages table.
                               [default: lineages.h5; required]

  -L, --length-buffer INTEGER  Fudge factor: a change in length less than this
                               much is not considered to reflect a cell
                               division.  [default: 5; required]

  -T, --trench-length TEXT     Length of a trench. Used to check whether a
                               cell is at the trench exit.  [required]

  --help                       Show this message and exit.
```

It is _highly_ recommended to correct trench numbering before performing lineage tracking (which is why the default input file is "newsegs.h5").

### write-kymographs

```
Usage: trenchcoat kymographs [OPTIONS]

  Pre-render intensity image kymographs & write to disk.

Options:
  -i, --images-file TEXT       Input HDF5 file with intensity images. Also
                               required for converting masks, because it
                               contains essential metadata.  [default:
                               HDF5/data.h5; required]

  -m, --masks-file TEXT        Input HDF5 file with segmentation masks.
  -R, --regions-file TEXT      HDF5 file with regions coordinates table.
                               Assumes that stage drift has been accounted for
                               and that trenches have been renumbered
                               ("corrected_trench_label").  [default:
                               corrected_regions.h5; required]

  -o, --out-dir TEXT           Directory to write HDF5 file with kymographs.
                               [required]

  -n, --num-cpu INTEGER RANGE  Number of CPUs to use.  [default: all CPUs]
  --help                       Show this message and exit.
```

### browse-kymographs

```
Usage: trenchcoat.py browse-kymographs [OPTIONS]

  Use Napari to browse kymographs. The regions file must be post-processed
  for correcting stage drift, and for merging cell information with trench
  information. Corrections are unimplemented.

Options:
  -i, --images-file TEXT      Input HDF5 file with images.  [default: HDF5/;
                              required]

  -m, --masks-file TEXT       Input HDF5 file with masks.
  -r, --regions-file TEXT     Input HDF5 file with regions (trenches must be
                              re-labeled for left/right drift).  [required]

  -L, --lineages-file TEXT    Input HDF5 file with cell measurements, re-
                              labeled trenches and cell lineages.

  -S, --settings-file TEXT    TrenchCoat settings file (YAML).  [default:
                              trenchcoat_settings.yaml; required]

  -V, --viewer-params TEXT    YAML file with image viewer params & information
                              for displaying column data.

  -d, --data-table-file TEXT  Path to HDF5 file containing cell measurements
                              for computed image.

  -P, --precomputed BOOLEAN   Specify whether kymographs of intensity images &
                              masks were pre-computed.  [required]

  --help                      Show this message and exit.
```



### Unimplemented or work-in-progress features:

#### lineage tracking: splitting merged cells

If two neighboring cells are erroneously merged, then the lineage tracking is disrupted. This is a non-trivial problem.

#### different trench detection algorithms

Add ability to pass in one of several different algorithms (update: seems to work, but hasn't been tested extensively).

#### parameter exploration using a GUI

Work in progress. The API changed recently, and a lot of work remains to make the interface more user friendly and work with more algorithms.

#### stage drift correction for trench detection

Work in progress. Will try to use cross-correlation between images to either further correct the stage drift, or as an independent calculation.

#### trench-measurements

Instead of segmenting cells within trenches, just measure the properties of entire trenches. This could be particularly useful at low magnification, when distinguishing between neighboring cells is difficult.

## An example workflow

1. Copy all ND2 files into a directory (named ``ND2``). These files will be analyzed together, using the same parameters.

2. Write a trenchcoat_settings.yaml file. This file defines the [viewer settings](https://napari.org/docs/dev/api/napari.html?highlight=viewer#napari.view_layers.view_image) for each microscopy ``channel``.

It also defines optional settings for flatfield corrections (works well) and registration (work in progress!).

In the future, it may be possible to auto-generate such a file, but for now the standard workflow is to copy-paste from previous files and to manually tweak them as needed. In practice, the most useful changes are to add, remove, or re-name the channels, their colormaps, names, and contrast limits.

For example:

    ```
    napari_settings:
        Phase:
            rgb            : False
            multiscale     : False
            colormap       : "gray"
            contrast_limits:
                - 500
                - 30000
            gamma          : 1.0
            interpolation  : "nearest"
            name           : "Phase contrast"
            opacity        : 1.0
            blending       : "additive"
            visible        : True

        YFP:
            rgb            : False
            multiscale     : False
            colormap       : "yellow"
            contrast_limits:
            - 0
            - 1000
            gamma          : 1.0
            interpolation  : "nearest"
            name           : "mVenus"
            opacity        : 1.0
            blending       : "additive"
            visible        : True

        CFP:
            rgb            : False
            multiscale     : False
            colormap       : "cyan"
            contrast_limits:
            - 0
            - 500
            gamma          : 1.0
            interpolation  : "nearest"
            name           : "SCFP3"
            opacity        : 1.0
            blending       : "additive"
            visible        : True

    # Optional. These are used for applying corrections on-the-fly.
    corrections:
        corrections_file: "corrections.h5"
        # Which image to correct which?
        # Input channel of image, return channel of correction
        flatfield_channel_mapping:
            Phase : NULL
            YFP   : "YFP"
            CFP   : "CFP"

        # Which image to correct which?
        # Input channel of image, return channel of correction
        registration_channel_mapping:
            Phase : NULL
            YFP   : "YFP"
            CFP   : "CFP"
	```

3. Inspect the ND2 files:

	```
	trenchcoat browse-nd2 -i ND2 -S trenchcoat_settings.yaml
	```

4. Convert the ND2 files to HDF5:

	```
	trenchcoat convert -i ND2 -o HDF5
	```

5. Confirm that the conversion was successful:

	```
	trenchcoat browse-hdf5 -i HDF5/data.h5 -S trenchcoat_settings.yaml
	```

6. Create a trench detection parameters file:

    ```
    algorithm: "algo_1"
    channel: "Phase-Fluor"
    # NOTE: Anything below "parameters" will be checked for its dimension, to enable auto-conversion from microns into pixels.
    parameters:
        trenches:
            tr_height:
                value: 25.0
                dimension: "microns"
            tr_width:
                value: 17
                dimension: "pixels"
            tr_spacing:
                value: 15.20
                dimension: "pixels"
            pad_left:
                value: 0
                dimension: "pixels"
            pad_right:
                value: 0
                dimension: "pixels"
            pad_top:
                value: 0
                dimension: "pixels"
            pad_bottom:
                value: 0
                dimension: "pixels"
    ```

7. Detect trenches and write the trench coordinates to disk:

	```
	trenchcoat trench-detect -i HDF5/data.h5 -o REGIONS -P trench_params.yaml --share-regions
	```

8. Create a segmentation parameters file:

	```
mKate2:
    algorithm: "fluor_sharpen"

    fluorescent_channel: MCHERRY"

    parameters:
        # -0.75 was pretty good, but occasionally cells seemed to be split up
        # about 1 frame too soon, and would merge again 1 frame.
        niblack_k: -0.75

        # Local window averaging
        # [9, 31] works well
        # First dim is along trench width, second is along trench length
        niblack_w: [13, 21]
        otsu_multiplier: 1.1

         # Important for avoiding background pixels in empty trenches
        garbage_otsu_value: 2200

        # Sharpening helps a lot!
        unsharp_mask_radius: 4
        unsharp_mask_amount: 5

        # Sometimes a small piece of a cell doesn't make it with the rest
        # Get rid of the fragments, lest they be considered separate cells
        # and interfere with the lineage tracking etc.
        min_size: 20.00
        # TODO Convert using px_mu? 5.00 # 60 pixels total; roughly 3.12 microns x 3.12 microns

	```

9. Segment cells and write their masks & measurements to disk:

	```
	trenchcoat segment -i HDF5/data.h5 -o SEGMENTATION -P seg_params.yaml -R REGIONS/regions_tables.h5
	```

10. View the microscopy images, trenches, and segmented cells:

	```
	trenchcoat browse-hdf5 -i HDF5/data.h5 -S napari_settings.yaml -m SEGMENTATION/MASKS/masks.h5 -r regions_tables.h5
	```

11. Generate kymographs & write to disk:
    ```
	trenchcoat write-kymographs -i HDF5/data.h5 -m segs.h5 -R corrected_regions.h5 -o KYMO/
	```

12. View kymographs (**NOTE**: double check this!):
    ```
    trenchcoat browse-kymographs -i KYMO_IMG/ -m KYMO_MASK/MASKS/masks.h5 -r corrected_regions.h5 -L lineages.h5 -S trenchcoat_settings.yaml
    ```
### Scripting

It is trivial to run a series of commands in a BASH script:

run.sh:

```
#!/bin/bash

trenchcoat convert -i ND2 -o HDF5

trenchcoat trench-detect -i HDF5/data.h5 -o REGIONS -P trench_params.yaml --share-regions

trenchcoat segment -i HDF5/data.h5 -o SEGMENTATION -P seg_params.yaml -R REGIONS/regions_tables.h5
```

In BASH:

```
# Make the file executable
chmod +x run.sh

# Run the analysis
./run.sh

# View the results
trenchcoat browse-hdf5 -i HDF5/data.h5 -S napari_settings.yaml -m SEGMENTATION/MASKS/masks.h5 -r REGIONS/regions_tables.h5
```

The use of BASH scripts also serves to document commands for future reference.

Likewise, it is possible to run many commands in SLURM, the HMS compute cluster:

```
#!/bin/bash

#SBATCH -c 8
#SBATCH -t 0-2:00
#SBATCH -p short
#SBATCH --mem=5000
#SBATCH -o hostname_%j.out
# Options are NONE, BEGIN, END, FAIL, REQUEUE, ALL, TIME_LIMIT_50, ARRAY_TASKS
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=example@hms.harvard.edu

###
#
# Things to keep in mind:
# - Don't forget to activate the TrenchCoat conda environment before sbatching the script to SLURM!
# - Make sure that my PATH variable contains ~/bin/trenchcoat.
# - Try running trenchcoat --help (maybe from an interactive session) to make sure it runs without errors.
#
# Arguments:
# 1. lane_number
# 2. trench_params_file_name
#
###

# For some reason, it can't find trenchcoat.
# Can we help with the PATH here? Even though it's clearly set when I start an interactive session...
PATH=$PATH:~/bin

nd2_dir="ND2"
hdf_dir="HDF5"
trenches_dir="REGIONS"
seg_dir="SEGMENTATION"

# NOTE remember to set this also to the same value for SBATCH -c
# NOTE can we leave this out, and have it automatically use all available CPUs?
num_cpu=8

# trench_params_1.3.yaml or trench_params_1.5.yaml
tr_params_file="../trench_params_1.3.yaml"

###

# Convert from ND2 to HDF5
trenchcoat convert -i "$nd2_dir" -o "$hdf_dir" -n $num_cpu

# Detect trenches
# NOTE odd-numbered lanes are 1.3 micron, even-numbered are 1.5 micron
trenchcoat trench-detect -i "$hdf_dir" -o "$trenches_dir" -n "$num_cpu" -P "$tr_params_file"

# Segment cells
trenchcoat segment -i "$hdf_dir" -o "$seg_dir" -P "../seg_params_trenches.yaml" -R "$trenches_dir/regions_tables.h5" -n "$num_cpu"

# Adjust trench numbering
trenchcoat correct-trench-numbering -i "$hdf_dir" -R "$trenches_dir" -M "$seg_dir" -O "relabeled_regions.h5" -m "relabeled_segmentations.h5" -H "$drift_method" -T $trench_threshold
```

## How does TrenchCoat work?

### ND2 to HDF5 conversion

A directory of ND2 files can be converted to the [HDF5 format](https://www.hdfgroup.org/solutions/hdf5/). This file format is not tied to Python, though TrenchCoat relies exclusively on the [Pytables](http://www.pytables.org) library.

#### Images

The intent is to convert once, and to *never* modify the HDF5 file. All subsequent operations (trench detection, cell segmentation, etc.) write to new files. It is good practice to set the HDF5 directory and its contents to read-only after converting.

- Image data are written as "chunked" arrays with lossless compression, resulting in roughly 50% reduction in file size.

- Images are *not* modified:
	+ No cropping is performed.
	+ No corrections are applied.

- Images are stored in a hierarchy:
	- ND2 file
		- FOV
			- Frame
				- Z (each Z-stack is stored as its own HDF5 file on disk)
					- Channel

- The decision to store images in HDF5 files by Z-stack was a compromise between having many small files on disk (bad) and allowing for parallel writes (good) across images which are not associated.

- All the HDF5 files are linked together into a "parent" HDF5 file. This strategy allows for independently writing to the Z-stack files (*via* standard filesystem paths), while still being able to afterwards browse the collection of files from a single HDF5 file (HDF5 does not support parallel writes). This strategy also allows for opening or transferring a single Z-stack HDF5 file independently, which could be useful for *e.g.* viewing in ImageJ.

##### A note about C arrays and Fortran arrays

There are two conventions for array indexing: C-style and Fortran-style arrays (see: [Row- and column-major order](https://en.wikipedia.org/wiki/Row-_and_column-major_order)). C-style arrays are more intuitive for relating array indices with memory positions and are the default in C & C++ (including [OpenCV](https://docs.opencv.org/2.4/index.html)), whereas Fortran-style arrays are more commonly used in mathematically inclined languages (*e.g.* Fortran, MATLAB, R, Julia, etc.). By default, Numpy uses C-style (row-major) indexing, but it also supports Fortran-style (column-major) indexing.

I decided to use column-major indexing for 2 reasons:

1. Interoperability with external libraries, in particular [ArrayFire](http://arrayfire.org/docs/index.htm) (for GPU acceleration; still a work in progress!)
2. Due to a quirk in the native ND2 format (see section about the modified nd2reader Python library)

In practice, the only effects that this decision has are that:

- images must be transposed when loading into Napari
- algorithms must use proper indexing

However, this decision might be unpopular for those wanting to share algorithms written using row-major indexing.

Options for dealing with this issue include:

- allowing the user to specify which indexing to use
- always using column-major indexing, and performing conversions as necessary
- writing a wrapper for functions to go between the two

My preference is to enforce a single indexing method so as to avoid the performance overhead of converting between the two.

#### Metadata

All metadata are copied and stored, with separate metadata for each ND2 file.

- Basic metadata, such as channel names and FOVs

- "Raw" metadata, containing information about acquisition settings of the microscope

- A table matching images to their acquisition times

### Trench detection (optional)

Trench detection is performed by analyzing phase contrast images. It is assumed that no rotational corrections need to be applied, and that trenches are vertically oriented.

After optionally cropping the image, an attempt is made to detect rows of trenches. (**NOTE:** row detection is still a work in progress! If it doesn't work well, then one option is to manually crop around a single row, and to repeat the analysis for each row with different cropping regions.)

Then, for each row, individual trenches are detected. The resulting coordinates are stored on disk in separate variable-length arrays.

Finally, all detected trench coordinates are copied to a global table, indexed by File, FOV, Frame, Z, and trench row number. The information in the variable-length arrays and the table are identical, but querying the table is much easier (and probably quicker) than traversing the arrays one at a time.

Trenches can either be detected for every time frame, or shared across time frames. Sharing trenches across frames eliminates the need to analyze every frame, but only works if there is little or no drift between frames.

Detection parameters are read from a YAML file. Currently there is no explicit support for alternative detection algorithms to be passed in as arguments, though this feature could be implemented in a manner similar to that used for cell segmentation (see below).

### Cell segmentation

A major goal for TrenchCoat is the support for arbitrary segmentation algorithms, with or without regions (*i.e.* trenches). To this end, a *generic workflow* was developed for segmenting Z-stacks.

#### Outline

The user must define:

- a *standardized, generic* segmentation algorithm wrapper function, which takes as input:
	+ the image stack
	+ the image stack -> index hashmap
	+ all parameters
- a *custom* segmentation algorithm function, which:
	+ takes any desired input (presumably these would be 1 or more images, and parameters)
	+ must return a segmentation mask, which may be a binary mask or a labeled mask, and must support 1 or more regions (*i.e.* 1 region: whole image. 2+ regions: trenches)

The results of each segmentation are stored on-disk in a new HDF5 file.

#### Details

##### User-defined settings

The user must specify, in a YAML file:

- how many segmentations to perform (each with a unique name)
- the segmentation algorithm (matching an algorithm defined above)
- which channels to use for segmentation (can be more than one)
- all parameters for that algorithm

##### Iteration

We iterate through every File, FOV, Frame, and Z-stack independently (but possibly in parallel).

NOTE: alternatively, a similar approach could be used to segment Frame stacks (kymographs), though this would require changing how iterations are performed and how image data are stored in memory.

##### Memory storage

A single array is used to store all relevant image data. The dimensions are (in column-major order) x - y - region - channel.

If regions (trenches) are provided, then the Z-stack is loaded while splitting the pixel data on a per-trench basis into a 4-D array. Otherwise, a 4-D array is still used, but the entire images are loaded into the zeroth position in the region dimension.

Note that this approach requires all regions (trenches) to have *identical* dimensions. One potential benefit of this approach is that applying the same operation across all regions can be further parallelized (*e.g.* on a GPU). One drawback is that there may be situations where trenches do not have identical dimensions. An alternative method would have to be written to take into account this possibility, where whole images are loaded, and then variable-dimensioned slices are looped & processed one at a time.

Each channel dimension is accompanied by a hashmap (dictionary) which relates its index in the array to its channel name. In this way, it is always possible to request a specific channel's image data. Thus, an image stack and its channel hashmap store the data in a manner accessible using string-based channel names and integer indices.

##### Data type

All images are loaded as 32-bit floating point (f32) images. This datatype was chosen because:

- f32 can store 16-bit unsigned integers (the default for our cameras) without loss of precision
- f32 supports operations, such as flat-field corrections, which apply non-integer changes
- f32 has enough precision for our needs (f64 is unnecessary)

##### Algorithms

When TrenchCoat first executes, all possible algorithms are defined in advance and stored in a hashmap (dictionary). In this way, any algorithm can be defined and referenced by a unique name. This name can be specified by the user in the YAML file (see above).

#### Example code for a segmentation algorithm

algorithms.py:

```
def run_single_threshold(stack, ch_to_index, params):
    """
    Very simple max thresholding example.
    """
    data = stack[ch_to_index[params["channel"]]]
    return single_threshold(data, **params["parameters"])

def single_threshold(data, cutoff):
    return data < cutoff
```

seg_params.yaml:

```
mSCFP3:
	algorithm: "single_threshold"
	channel: "CFP"

	parameters:
		cutoff: 500
```

Or, a more complex segmentation algorithm:

algo_sharp.py:

```
def run_fluor_sharp_segmentation(stack, ch_to_index, params):
    """
    A version which uses thresholding on fluor & phase channels
    Thresholding which also uses Niblack (mean, std dev) & watershedding
    """
    stack_fl = stack[..., ch_to_index[params["fluorescent_channel"]]]

    return fluor_sharpen_segmentation(stack_fl, **params["parameters"])

def fluor_sharpen_segmentation(
    stack_fl,
    otsu_multiplier,
    garbage_otsu_value,
    niblack_w,
    niblack_k,
    min_size,
    unsharp_mask_radius,
    unsharp_mask_amount
):
    """
    Use fluorescence.
    Key difference from other Niblack-based routines is the use of sharpening.
    """
    # Sharpen the fluorescence image
    for z in range(stack_fl.shape[2]):
        stack_fl[..., z] = unsharp_mask(
            stack_fl[..., z],
            radius=unsharp_mask_radius,
            amount=unsharp_mask_amount,
            preserve_range=True
        )

    # Threshold is calculated using Otsu method, which samples the distribution
    # of pixel intensities, and then scaled by the otsu multiplier, which might
    # differ depending on empirical experience with a given channel or dataset.
    #
    # Then set pixels < threshold to zero, but preserve values > threshold.
    # Helps with the mean + stdev*k calculations in Niblack.
    # It is "bimodal" because the values are either zero, or they are the
    # original value.
    # NOTE: is the broadcasting operating as expected for threshold[z] vs bimodal[..., z]?
    # bimodal = ((stack_fl >= threshold) * stack_fl).astype(stack_fl.dtype, casting="safe", copy=False)
    bimodal = numpy.empty(stack_fl.shape, dtype=numpy.float32, order="F")
    for z in range(bimodal.shape[2]):
        threshold = threshold_otsu(stack_fl[..., z]) * otsu_multiplier

        bimodal[..., z] = ((stack_fl[..., z] >= threshold) \
                          * stack_fl[..., z]).astype(stack_fl.dtype, casting="safe", copy=False)

        # Without using an if statement, multiply: bimodal *= (threshold > garbage_otsu_value)
        # Either x 1, stays the same, or x 0, and then all subsequent steps do nothing.
        bimodal[..., z] *= threshold > garbage_otsu_value

    # Apply Niblack method
    niblack = numpy.empty(stack_fl.shape, dtype=numpy.float32, order="F")
    for z in range(niblack.shape[2]):
        niblack[..., z] = threshold_niblack(
                            bimodal[..., z],
                            window_size=niblack_w,
                            k=niblack_k
                        )

    # Only keep pixels which are < than the mean*stdev*k (in rectangular region, size w).
    # NOTE does it help to specify it explicitly, with F ordering?
    # or does it actually make things worse, because the next line then trashes it and makes a new array???
    bimodal = bimodal.astype(numpy.bool, casting="unsafe", copy=False)
    mask = (stack_fl > niblack) * bimodal

    markers = numpy.empty(stack_fl.shape, dtype=numpy.uint16, order="F")
    for z in range(mask.shape[2]):
        # Fill in internal holes in the "mask"
        mask[..., z] = binary_fill_holes(mask[..., z])

        # Remove small objects.
        # Connectivity of 1 seems to help a lot in preventing adjacent cells from being merged?
        # But connectivity of 2 can help if the eroded bits are very small.
        mask[..., z] = remove_small_objects(
            mask[..., z],
            min_size=min_size,
            connectivity=1
        )

        # Label the regions
        # NOTE is there a way to change the labeling so that it
        # coincides with the "top down" view that we associate with trenches?
        # (it seems to go more left-to-right, which doesn't make sense!)
        # Yes, by making an intermediate, tranposed labels array, and then
        # transposing it back & storing it again. Could maybe improve
        # performance by tweaking their algo. to operate in the correct dimensional order?
        markers[..., z] = label(
            mask[..., z].T,
            connectivity=1
        ).T

    del niblack
    del bimodal


    # Make a basin (invert max, min pixels). Use the original intensity image for
    # the pixel values to use for the gradient ascent. Take advantage of the mask
    # to exclude undesired pixels.
    # NOTE it would be nice if the watershed algo. allowed for a max priority queue.
    # Would save us a step here.
    # NOTE it shouldn't matter if we use the max across all regions, or a same region,
    # it only matters that it's truly larger than all other values.
    image = stack_fl.max() - stack_fl

    # Store a stack of labeled segmentation masks
    # Must be integer type (what is returned by watershedding, i.e. integer labels.)
    result = numpy.empty(stack_fl.shape, dtype=numpy.uint16, order="F")
    for z in range(result.shape[2]):
        # FIXME this fixes really odd behavior where
        # the watershed refuses to handle the mask properly
        # unless I multiply it by a bunch of boolean ones.
        o = numpy.ones(mask[..., z].shape, dtype=numpy.bool)
        result[..., z] = watershed(
                            image=image[..., z],
                            markers=markers[..., z],
                            mask=mask[..., z] * o # multiply by boolean 1s
                         )

    return result
```

seg_params.yaml:

```
mKate2:
    algorithm: "fluor_sharpen"

    fluorescent_channel: MCHERRY"

    parameters:
        # -0.75 was pretty good, but occasionally cells seemed to be split up
        # about 1 frame too soon, and would merge again 1 frame.
        niblack_k: -0.75

        # Local window averaging
        # [9, 31] works well
        # First dim is along trench width, second is along trench length
        niblack_w: [13, 21]
        otsu_multiplier: 1.1

         # Important for avoiding background pixels in empty trenches
        garbage_otsu_value: 2200

        # Sharpening helps a lot!
        unsharp_mask_radius: 4
        unsharp_mask_amount: 5

        # Sometimes a small piece of a cell doesn't make it with the rest
        # Get rid of the fragments, lest they be considered separate cells
        # and interfere with the lineage tracking etc.
        min_size: 20.00
        # TODO Convert using px_mu? 5.00 # 60 pixels total; roughly 3.12 microns x 3.12 microns
```

Note that the name of the segmentation ("mSCFP3", "mKate2") needn't match the name of the channel ("CFP", "MCHERRY").

**The key feature here is that the function signatures are _always_ the same!! This allows us to pass in arbitrary segmentation functions!!**
```run_function_name(stack, ch_to_index, params):```


A brief explanation of these examples:
- The first function (run\_single\_threshold, or run\_fluor\_sharp\_segmentation) is a wrapper function. It takes as input the image stack, the dictionary mapping channel names to indices within the image stack, and a dictionary of parameters. All functions *must* share this function signature.

- This function then specifies which channel is the one to be used for segmentation by interrogating the "channel" parameter, and uses the ch\_to_index dictionary to determine which image to pass onwards.

- The second function (single\_threshold, fluor\_sharp\_segmentation) is then invoked using the image and parameters, and its result is returned back to the first function, which in turn returns it unchanged. In the first case (single\_threshold), a binary mask will be returned after determining which pixels are less than 500.

- The first example does not use loops to process individual regions, because comparison operators are automatically applied across arrays of arbitrary dimensionality The second example explicitly loops across array.shape\[2\], which is the dimension of regions (trenches).

More complex examples are to be found in ``algorithms.py``.

### <a id="nd2reader">The modified nd2reader library</a>

The original nd2library [contains](https://github.com/jimrybarski/nd2reader/blob/master/nd2reader/driver/v3.py) the following comment:

```
# The images for the various channels are interleaved within the same array. For example, the second image
# of a four image group will be composed of bytes 2, 6, 10, etc. If you understand why someone would design
# a data structure that way, please send the author of this library a message.
```

As a result, whenever a request is made to read a particular image of a given channel, all channels are loaded, and then all but one of the channels *are discarded*. For an image with 4 channels, iterating across all 4 channels requires reading the same data from disk 4 times! This problem is amplified in proportion to the number of channels.

To resolve this problem, I made added a function which allows requesting the entire stack of channels in a single pass.

The following two additions:

main.py:

```
def get_image_stack(self, frame_number, field_of_view, z_level):
        return self._parser.driver.get_image_stack(frame_number,
                                                   field_of_view,
                                                   z_level,
                                                   self.height,
                                                   self.width)
```

driver/v3.py:

```
def get_image_stack(self, frame_number, field_of_view, z_level, height, width):
    """
    Attempts to get Image stack based on attributes alone.
    :type frame_number:  int
    :type field_of_view: int
    :type z_level:       int
    :type height:        int
    :type width:         int
    :rtype: Fortran-ordered Numpy array, or None
    # NOTE: dict mapping channel name to index will have to be inferred from somewhere else!!
    """
    # Use the chunk map to find the group of images, and read into a flat sequence of intensities
    image_group_number = self._calculate_image_group_number(frame_number, field_of_view, z_level)
    chunk = self._label_map.get_image_data_location(image_group_number)
    data = read_chunk(self._file_handle, chunk)

    # Convert the flat intensities to an "array"
    # NOTE: why not an ndarray?
    image_group_data = array.array("H", data)

    # NOTE What are the first 4 bytes? Timestamp? Is it bytes, or words, or ... ?
    image_data_start = 4

    # Based on the total length of the data, and the known width and height,
    # extrapolate the number of channels in the data.
    # NOTE why not use integer divide, // ?
    number_of_true_channels = (len(image_group_data) - image_data_start) // (height * width)
    #number_of_true_channels = int((len(image_group_data) - image_data_start) / (height * width))

    # Zeroth dimension is the channel
    # zeroth index -> first channel, 1st index -> 2nd channel etc.
    # This produces an array with channels as the first of 3 indices, which is inefficient for F-ordering.
    reshaped = np.reshape(image_group_data[image_data_start:], (number_of_true_channels, width, height), order='F')

    # Fortran ordering works best when the last index is the slowest-changing
    # For image stacks, this typically means that 2D slices are stored in the 1st 2 dims,
    # and therefore a whole slice is in the 3rd dimension.
    # NOTE because the data are Fortran-ordered, they will appear Transposed
    # when viewing with matplotlib pyplot. However, it's more appropriate to Transpose later,
    # and to leave the data intact now.
    image_stack = np.empty((width, height, number_of_true_channels), dtype=reshaped.dtype, order='F')
    for c in range(number_of_true_channels):
        image_stack[..., c] = reshaped[c]


    # Skip images that are all zeros! This is important, since NIS Elements creates blank "gap" images if you
    # don't have the same number of images each cycle. We discovered this because we only took GFP images every
    # other cycle to reduce phototoxicity, but NIS Elements still allocated memory as if we were going to take
    # them every cycle.
    if np.any(image_stack):
        # NOTE: does this mean the first 8 bytes are the timestamp? 4 bytes? Why not parse out the timestamp first, then load the rest into the array? This was a goofy way of doing things by Rybarski...
        # NOTE do we care about the timestamp??
        #timestamp = struct.unpack("d", data[:8])[0]
        #return timestamp, image_stack
        return image_stack
    raise NoImageError
```

Observations:

- It turns out that, when using 'F' ordering, splitting the interleaved images by channel is straightforward. I'm guessing that, internally, ND2 files are somehow created using column-major ordering.
- I have not actually benchmarked to what extent this addition improves conversion speeds.
- This added function only applies when reading ND2 files, and has no effect when later working with HDF5 files.
