#!/bin/sh
OUTPUT_FOLDER=/n/groups/paulsson/paulsson-home/conda-channel
for package in guppy; do
	conda build --output-folder "$OUTPUT_FOLDER" "$package"
done
