#!/bin/bash

set -ex

# use out-of-tree build
mkdir -p _build
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

# handle cross compiling
if [[ $build_platform != $target_platform ]]; then
	CONFIGURE_ARGS="${CONFIGURE_ARGS} --disable-help2man"
fi

# only link libraries we actually use
export CFITSIO_LIBS="-L${PREFIX}/lib -lcfitsio"
export GSL_LIBS="-L${PREFIX}/lib -lgsl"

# configure
${SRC_DIR}/configure \
	--disable-doxygen \
	--enable-cfitsio \
	--prefix="${PREFIX}" \
	${CONFIGURE_ARGS} \
;

# build
make -j ${CPU_COUNT} V=1 VERBOSE=1

# test
if [[ "${build_platform}" == "${target_platform}" ]]; then
	make -j ${CPU_COUNT} V=1 VERBOSE=1 check
fi

# install
make -j ${CPU_COUNT} V=1 VERBOSE=1 install
