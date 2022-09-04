#!/bin/bash

set -ex

# use out-of-tree build
mkdir -pv _build
cd _build

# customisation for LALSuite development CI
if [[ "${GITLAB_CI}" == "true" ]] && [[ "x${CI_COMMIT_TAG}" == x ]]; then
	# declare nightly builds
	if [ "${CI_PIPELINE_SOURCE}" = "schedule" ] || [ "${CI_PIPELINE_SOURCE}" = "web" ]; then
		CONFIGURE_ARGS="${CONFIGURE_ARGS} --enable-nightly"
	fi
# production builds ignore GCC warnings
else
	CONFIGURE_ARGS="${CONFIGURE_ARGS} --disable-gcc-flags"
fi

# select FFT implementation
if [[ "${fft_impl}" == "mkl" ]]; then
	CONFIGURE_ARGS="${CONFIGURE_ARGS} --disable-static --enable-intelfft"
fi

# only link libraries we actually use
export GSL_LIBS="-L${PREFIX}/lib -lgsl"
export HDF5_LIBS="-L${PREFIX}/lib -lhdf5 -lhdf5_hl"

# configure
${SRC_DIR}/configure \
	--disable-doxygen \
	--disable-help2man \
	--disable-python \
	--disable-swig-octave \
	--disable-swig-python \
	--enable-swig-iface \
	--prefix="${PREFIX}" \
	${CONFIGURE_ARGS} \
;

# build
make -j ${CPU_COUNT} V=1 VERBOSE=1 HDF5_LIBS="${HDF5_LIBS}"

# test
if [[ "${build_platform}" == "${target_platform}" ]]; then
	make -j ${CPU_COUNT} V=1 VERBOSE=1 check
fi
