#!/bin/bash
#
# Configure, built, test, and install the Python language bindings
# for a LALSuite subpackage.
#

set -e
pushd ${SRC_DIR}

# only link libraries we actually use
export GSL_LIBS="-L${PREFIX}/lib -lgsl"
export CFITSIO_LIBS="-L${PREFIX}/lib -lcfitsio"

# configure only python bindings and pure-python extras
./configure \
	--prefix=$PREFIX \
	--disable-swig-iface \
	--enable-swig-python \
	--enable-python \
	--disable-doxygen \
	--disable-gcc-flags \
	--enable-silent-rules || { cat config.log; exit 1; }

# build
make -j ${CPU_COUNT} -C swig
make -j ${CPU_COUNT} -C python

# test
make -j ${CPU_COUNT} -C test check

# install
make -j ${CPU_COUNT} -C swig install-exec-am  # swig bindings
make -j ${CPU_COUNT} -C python install  # pure-python extras

popd
