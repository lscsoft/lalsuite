#!/bin/bash

set -ex

# enable nightly mode for CI
if [ "${CI_PIPELINE_SOURCE}" = "schedule" ] || [ "${CI_PIPELINE_SOURCE}" = "web" ]; then
	ENABLE_NIGHTLY="--enable-nightly"
fi

# when running on gitlab-ci, we are not using a production
# build, so we don't want to use NDEBUG
export CPPFLAGS="${CPPFLAGS} -UNDEBUG"

# only link libraries we actually use
export CFITSIO_LIBS="-L${PREFIX}/lib -lcfitsio"
export GSL_LIBS="-L${PREFIX}/lib -lgsl"

# configure
./configure \
	--disable-doxygen \
	--enable-cfitsio \
	--enable-help2man \
	--prefix=${PREFIX} \
	${ENABLE_NIGHTLY} \
;

# build
make -j ${CPU_COUNT} V=1 VERBOSE=1

# test
make -j ${CPU_COUNT} V=1 VERBOSE=1 check

# install
make -j ${CPU_COUNT} V=1 VERBOSE=1 install
