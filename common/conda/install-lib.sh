#!/bin/bash

set -ex

_make="make -j ${CPU_COUNT} V=1 VERBOSE=1"

cd _build

# install library and headers
${_make} -C lib install

# install SWIG binding definitions and headers
${_make} -C swig install-data

# install pkg-config
${_make} install-pkgconfigDATA
