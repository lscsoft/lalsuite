#!/bin/bash

set -e

# select FFT implementation
if [[ "${fft_impl}" == "mkl" ]]; then
    FFT_CONFIG_ARGS="--disable-static --enable-intelfft"
else
    FFT_CONFIG_ARGS=""
fi

# only link libraries we actually use
export GSL_LIBS="-L${PREFIX}/lib -lgsl"

./configure \
	--prefix="${PREFIX}" \
	--enable-swig-iface \
	--disable-swig-octave \
	--disable-swig-python \
	--disable-python \
	--disable-gcc-flags \
	--enable-silent-rules \
	${FFT_CONFIG_ARGS}
make -j ${CPU_COUNT}
make -j ${CPU_COUNT} check
make install
