# snakenbake

## Installation
After cloning this repository, `cd snakenbake` and run the following (replace `snakenbake` with the desired name for your conda environment):
```
conda env create -n snakenbake -f environment.yml
echo "conda activate snakenbake" > .envrc
direnv allow
curl http://ftp.gnome.org/pub/GNOME/sources/ttf-bitstream-vera/1.10/ttf-bitstream-vera-1.10.tar.bz2|tar xvj
```

The last command downloads the fonts, which are required for text rendering.