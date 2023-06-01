#!/bin/bash

# Run a series of Python programs to analyze a mother machine experiment

# 1. Convert the ND2 file to HDF5

ND2_file_path="/home/andrew/Martens/2019-12-04 Delay Time/delaytime_9_16_100x_001.nd2"
HDF5_output_dir="/home/andrew/Martens/2019-12-04 Delay Time/HDF5_0_29"

# Using 8 processes gave a 4.24x speedup over using 1 process.
# Transferring bandwidth for 4 time frames at once would require ~1 Gbps alone!
# Not surprising that this saturates quickly.
# Bandwidth limitation? 4 is probably a reasonable value.
num_cpu=4
min_frame_number=0
max_frame_number=29 # 163

# NOTE that it is essential to put quotes around the path names here, in case there are spaces!
./trenchcoat.py convert -i "$ND2_file_path" -o "$HDF5_output_dir" -n $num_cpu -m $min_frame_number -M $max_frame_number

# 2. Detect trenches
# TODO first, have it detect trench rows as a separate step? -> to get crop_top & crop_bottom values
# (possibly multiple such values). And these values might change by FOV?
# TODO to mitigate the possibility of HDF5 corruption & having to re-convert the ND2 file,
# move to putting all tables, masks etc. in separate files & link them into the FOV files.

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

./trenchcoat.py trench-detect -i "$HDF5_dir" -p "$parent_name" -n $num_cpu -d $min_distance -c $cutoff -W $half_width -T $crop_top -B $crop_bottom -I "$identifier" -C "$channel"

# 3. Analyze whole trenches

HDF5_dir="/home/andrew/Martens/2019-12-04 Delay Time/HDF5_0_29/"
parent_name="delaytime_9_16_100x_001.h5"
num_cpu=1
identifier="1"
crop_top=1000
crop_bottom=1325

./trenchcoat.py trench-measurements -i "$HDF5_dir" -p "$parent_name" -n $num_cpu -I "$identifier" -T $crop_top -B $crop_bottom

# # 4. Define bounding boxes within each trench
# # NOTE: this was less useful that I had hoped (help provide constraints for cell segmentation), and is error-prone.
# # Segmentation alone works well enough.
#
# HDF5_dir="/home/andrew/Martens/2019-12-04 Delay Time/HDF5_0_29/"
# parent_name="delaytime_9_16_100x_001.h5"
# num_cpu=1
# identifier="1"
# crop_top=1000
# crop_bottom=1325
# min_peak_distance=20
# prominence_file="prominence_values.tsv"
#
# ./cell_detection.py -i "$HDF5_dir" -p "$parent_name" -n $num_cpu -I "$identifier" -T $crop_top -B $crop_bottom -d $min_peak_distance -P $prominence_file

# 5. Segment cells within each trench

HDF5_dir="/home/andrew/Martens/2019-12-04 Delay Time/HDF5_0_29/"
parent_name="delaytime_9_16_100x_001.h5"
num_cpu=8
identifier="1"
crop_top=1000
crop_bottom=1325
params_file="params.txt"

./trenchcoat.py segment -i "$HDF5_dir" -p "$parent_name" -n $num_cpu -I "$identifier" -T $crop_top -B $crop_bottom -P $params_file

# 6. Generate kymographs
in_dir="/home/andrew/Martens/2019-12-04 Delay Time/HDF5_0_29/"
out_file="kymographs"
num_cpu=8
identifier="1"
crop_top=1000
crop_bottom=1325

./trenchcoat write-kymographs -i "$in_dir" -o "$out_file" -n $num_cpu -I $identifier -T $crop_top -B $crop_bottom
