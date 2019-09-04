#!/bin/bash

set -e

# only link libraries we actually use
export GSL_LIBS="-L${PREFIX}/lib -lgsl"
export CFITSIO_LIBS="-L${PREFIX}/lib -lcfitsio"

./configure \
	--prefix="${PREFIX}" \
	--enable-swig-iface \
	--disable-swig-octave \
	--disable-swig-python \
	--disable-python \
	--disable-gcc-flags \
	--enable-silent-rules \
	--enable-cfitsio
make -j ${CPU_COUNT}
make -j ${CPU_COUNT} check
make install
