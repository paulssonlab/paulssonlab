mkdir $PREFIX/bin
mkdir $PREFIX/lib
mkdir $PREFIX/data
cp bin/guppy_* $PREFIX/bin
cp bin/minimap2 $PREFIX/bin
cp lib/*.so* $PREFIX/lib
cp -R data/* $PREFIX/data
