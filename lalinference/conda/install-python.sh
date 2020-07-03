#!/bin/bash
#
# Configure, built, test, and install the Python language bindings
# for a LALSuite subpackage.
#

set -e

# when running on gitlab-ci, we are not using a production
# build, so we don't want to use NDEBUG
export CPPFLAGS="${CPPFLAGS} -UNDEBUG"

# configure only python bindings and pure-python extras
./configure \
	--disable-doxygen \
	--disable-swig-iface \
	--enable-help2man \
	--enable-mpi \
	--enable-python \
	--enable-swig-python \
	--prefix=$PREFIX \
;

# swig bindings
make -j ${CPU_COUNT} V=1 VERBOSE=1 -C swig
make -j ${CPU_COUNT} V=1 VERBOSE=1 -C swig install-exec-am

# python modules
make -j ${CPU_COUNT} V=1 VERBOSE=1 -C python
make -j ${CPU_COUNT} V=1 VERBOSE=1 -C python install

# python scripts
make -j ${CPU_COUNT} V=1 VERBOSE=1 -C bin bin_PROGRAMS="" dist_bin_SCRIPTS=""
make -j ${CPU_COUNT} V=1 VERBOSE=1 -C bin bin_PROGRAMS="" dist_bin_SCRIPTS="" install
