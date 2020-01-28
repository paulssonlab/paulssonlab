#!/bin/bash

# Run a series of Python programs to analyze a mother machine experiment

# 1. Convert the ND2 file to HDF5

# TODO: determine whether running in WSL if faster or slower than running directly on Windows.
# Under WSL, takes ~0.13 sec / time frame when using 4 processes
# (8344 sec for the entire 100x163 file: 163 time frames, 100 FOV, 4 images / time frame = 65,200 total images to convert)

ND2_file_path="/home/andrew/Martens/2019-12-04 Delay Time/delaytime_9_16_100x_001.nd2"
HDF5_output_dir="/home/andrew/Martens/2019-12-04 Delay Time/HDF5_0_29"

# Using 8 processes gave a 4.24x speedup over using 1 process.
# Transferring bandwidth for 4 time frames at once would require ~1 Gbps alone!
# Not surprising that this saturates quickly.
# Bandwidth limitation? 4 is probably a reasonable value.
num_cpu=4

# Note is this inclusive or not? What if the range is out of range?
# 0-2 -> Frame_0 and Frame_1
min_frame_number=0
max_frame_number=29 # 163

# NOTE that it is essential to put quotes around the path names here, in case there are spaces!
./convert.py -i "$ND2_file_path" -o "$HDF5_output_dir" -n $num_cpu -m $min_frame_number -M $max_frame_number

# 2. Detect trenches
# TODO first, have it detect trench rows as a separate step? -> to get crop_top & crop_bottom values
# (possibly multiple such values). And these values might change by FOV?

HDF5_dir="/home/andrew/Martens/2019-12-04 Delay Time/HDF5_0_29/"
parent_name="delaytime_9_16_100x_001.h5"
num_cpu=8
min_distance=25
cutoff=0.35
half_width=15
crop_top=1000
crop_bottom=1325
identifier="1"
channel="BF"

./trench_detect.py -i "$HDF5_dir" -p "$parent_name" -n $num_cpu -d $min_distance -c $cutoff -W $half_width -T $crop_top -B $crop_bottom -I "$identifier" -C "$channel"

# 3. Analyze whole trenches

HDF5_dir="/home/andrew/Martens/2019-12-04 Delay Time/HDF5_0_29/"
parent_name="delaytime_9_16_100x_001.h5"
num_cpu=8
identifier="1"
crop_top=1000
crop_bottom=1325

./trench_measurements.py -i "$HDF5_dir" -p "$parent_name" -n $num_cpu -I "$identifier" -T $crop_top -B $crop_bottom 

# 4. Segment cells within each trench

HDF5_dir="/home/andrew/Martens/2019-12-04 Delay Time/HDF5/"
parent_name="delaytime_9_16_100x_001.h5"
num_cpu=8 # How many can be run in parallel?
identifier="1"
crop_top=1000
crop_bottom=1325
params_file="params.txt"

./segmentation.py -i "$HDF5_dir" -p "$parent_name" -n $num_cpu -I "$identifier" -T $crop_top -B $crop_bottom -P $params_file

# 5. --> Look at data using a Jupyter notebook?

