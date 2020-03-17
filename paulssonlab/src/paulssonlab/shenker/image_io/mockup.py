### Principles

# each raw pixel has a physical address
# each resampled pixel has a physical address to the output of the resampling
#   resampling metadata can be used to map one transformed pixel address to (many) raw pixel physical addresses
# (non-resampled) transformations define a mapping from transformed address to physical address

# how does transformed address structure impact how physical pixel data is accessed?
#   what if I want trench crop cubes over time, but pixel data is stored in ND2?
#   how can I configure readahead buffering?
#   on local disk and in memory caching of both whole frames and of trench cube stacks?

# segmenting in cell crops or trench crops, need to reproject to whole frames

## Implementation 1:
# Physical pixel address: blob (exp ID)/t/ch/y,x
# Logical Experiment address: combine multiple ND2s into the same expID/t
# Join across trench positions (need to find chip coordinate system)
# need to inform caching policy based on accesses at the logical level

### List

# single ND2

# folder of TIFFs

# folder of ND2s

# zarrs (+different ways to chunk, divide hierarchy between zarr groups and ndim arrays)

### Image/metadata caching (caching in memory, on local disk, remote)

### Copy/Conversions

# copy of zarr -> zarr (but mapping fov/t/channel/y,x -> fov/channel/trench/t,y,x)

# watch single ND2 -> copy to zarr (SSH/GCP)

# watch folder for TIFFs -> copy to zarr (SSH/GCP)

# watch folder for ND2s -> copy to zarr (SSH/GCP)

# watch single ND2 -> crop trenches -> copy to zarr (SSH/GCP)

# watch folder of ND2s -> copy to zarr (SSH/GCP) -> crop trenches

### Image cropping/transformations

# image -> crop trenches -> cell crops

# two grids/experiments/microscopes (same chip) -> crop trenches -> cell crops with masks
# image -> derotate (?)

### How to store trench metadata
### how to store cell metadata
### how to store tabular data (and tie to physical pixels, or transformed pixels)
### how to store metadata after a processing step
### how to store lineage info

### Analysis
# trench crops over time

# cell segments

# cell segments over lineage history
