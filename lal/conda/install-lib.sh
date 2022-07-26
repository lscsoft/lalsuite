#!/bin/bash

set -ex

cd _build

# install library and headers
make -j ${CPU_COUNT} V=1 VERBOSE=1 -C lib install HDF5_LIBS="-L${PREFIX} -lhdf5 -lhdf5_hl"

# install SWIG binding definitions and headers
make -j ${CPU_COUNT} V=1 VERBOSE=1 -C swig install-data

# install pkg-config
make -j ${CPU_COUNT} V=1 VERBOSE=1 install-pkgconfigDATA
