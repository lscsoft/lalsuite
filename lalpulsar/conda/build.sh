#!/bin/bash

set -e

# when running on gitlab-ci, we are not using a production
# build, so we don't want to use NDEBUG
export CPPFLAGS="${CPPFLAGS} -UNDEBUG"

# only link libraries we actually use
export GSL_LIBS="-L${PREFIX}/lib -lgsl"
export CFITSIO_LIBS="-L${PREFIX}/lib -lcfitsio"

./configure \
	--prefix="${PREFIX}" \
	--enable-swig-iface \
	--disable-swig-octave \
	--disable-swig-python \
	--disable-python \
	--enable-silent-rules \
	--enable-cfitsio
make -j ${CPU_COUNT}
make -j ${CPU_COUNT} check
make install
