#!/bin/bash

set -ex

# use out-of-tree build
mkdir -pv _build
cd _build

# enable nightly mode for CI
if [ "${CI_PIPELINE_SOURCE}" = "schedule" ] || [ "${CI_PIPELINE_SOURCE}" = "web" ]; then
	ENABLE_NIGHTLY="--enable-nightly"
fi

# when running on gitlab-ci, we are not using a production
# build, so we don't want to use NDEBUG
export CPPFLAGS="${CPPFLAGS} -UNDEBUG"

# select FFT implementation
if [[ "${fft_impl}" == "mkl" ]]; then
    FFT_CONFIG_ARGS="--disable-static --enable-intelfft"
else
    FFT_CONFIG_ARGS=""
fi

# only link libraries we actually use
export GSL_LIBS="-L${PREFIX}/lib -lgsl"
export HDF5_LIBS="-L${PREFIX}/lib -lhdf5 -lhdf5_hl"

# configure
${SRC_DIR}/configure \
	--disable-doxygen \
	--disable-python \
	--disable-swig-octave \
	--disable-swig-python \
	--enable-help2man \
	--enable-swig-iface \
	--prefix="${PREFIX}" \
	${FFT_CONFIG_ARGS} \
	${ENABLE_NIGHTLY} \
;

# build
make -j ${CPU_COUNT} V=1 VERBOSE=1 HDF5_LIBS="${HDF5_LIBS}"

# test
make -j ${CPU_COUNT} V=1 VERBOSE=1 check
