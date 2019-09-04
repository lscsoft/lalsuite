#!/bin/bash
#
# Configure, built, test, and install the Python language bindings
# for a LALSuite subpackage.
#

set -e
pushd ${SRC_DIR}

# if we're using MKL, the C library will have been built with
# --enable-intelfft, so we have to use that here as well
if [[ "${fft_impl}" == "mkl" ]]; then
    FFT_CONFIG_ARGS="--disable-static --enable-intelfft"
else
    FFT_CONFIG_ARGS=""
fi

# only link libraries we actually use
export GSL_LIBS="-L${PREFIX}/lib -lgsl"

# configure only python bindings and pure-python extras
./configure \
	--prefix=$PREFIX \
	--disable-swig-iface \
	--enable-swig-python \
	--enable-python \
	--disable-doxygen \
	--disable-gcc-flags \
	--enable-silent-rules \
	${FFT_CONFIG_ARGS} \
|| { cat config.log; exit 1; }

# build
make -j ${CPU_COUNT} -C swig
make -j ${CPU_COUNT} -C python

# install
make -j ${CPU_COUNT} -C swig install-exec-am  # swig bindings
make -j ${CPU_COUNT} -C python install  # pure-python extras

popd
