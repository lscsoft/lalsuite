#!/bin/bash

set -ex

# use out-of-tree build
mkdir -pv _build
cd _build

# customisation for LALSuite development CI
if [[ "${GITLAB_CI}" == "true" ]] && [[ -z "${CI_COMMIT_TAG+x}" ]]; then
	# allow debugging information
	export CPPFLAGS="${CPPFLAGS} -UNDEBUG"

	# declare nightly builds
	if [ "${CI_PIPELINE_SOURCE}" = "schedule" ] || [ "${CI_PIPELINE_SOURCE}" = "web" ]; then
		CONFIGURE_ARGS="${CONFIGURE_ARGS} --enable-nightly"
	fi
# production builds ignore GCC warnings
else
	CONFIGURE_ARGS="${CONFIGURE_ARGS} --disable-gcc-flags"
fi

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
make -j ${CPU_COUNT} V=1 VERBOSE=1

# test
if [[ "${build_platform}" == "${target_platform}" ]]; then
	make -j ${CPU_COUNT} V=1 VERBOSE=1 check
fi
