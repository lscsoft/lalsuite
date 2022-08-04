#!/bin/bash
#
# Configure, built, test, and install the Python language bindings
# for a LALSuite subpackage.
#

set -ex

# build python in a sub-directory using a copy of the C build
_builddir="_build${PY_VER}"
cp -r _build ${_builddir}
cd ${_builddir}

# customisation for LALSuite development CI
if [[ "${GITLAB_CI}" == "true" ]] && [[ "x${CI_COMMIT_TAG}" == x ]]; then
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
if [[ "${build_platform}" != "${target_platform}" ]]; then
	# help2man doesn't work when cross compiling
	CONFIGURE_ARGS="${CONFIGURE_ARGS} --disable-help2man"
fi

# if we're using MKL, the C library will have been built with
# --enable-intelfft, so we have to use that here as well
if [[ "${fft_impl}" == "mkl" ]]; then
    CONFIGURE_ARGS="${CONFIGURE_ARGS} --disable-static --enable-intelfft"
fi

# only link libraries we actually use
export GSL_LIBS="-L${PREFIX}/lib -lgsl"

# configure only python bindings and pure-python extras
${SRC_DIR}/configure \
	--disable-doxygen \
	--disable-swig-iface \
	--enable-python \
	--enable-swig-python \
	--prefix="${PREFIX}" \
	${CONFIGURE_ARGS} \
;

# build
make -j ${CPU_COUNT} V=1 VERBOSE=1 -C swig
make -j ${CPU_COUNT} V=1 VERBOSE=1 -C python

# test
if [[ "${build_platform}" == "${target_platform}" ]]; then
	make -j ${CPU_COUNT} V=1 VERBOSE=1 check -C swig
fi

# install
make -j ${CPU_COUNT} V=1 VERBOSE=1 -C swig install-exec  # swig bindings
make -j ${CPU_COUNT} V=1 VERBOSE=1 -C python install  # pure-python extras
