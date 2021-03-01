#!/bin/bash

set -ex

# install from python build directory
_pybuilddir="_build${PY_VER}"
cd ${_pybuilddir}

# install binaries
make -j ${CPU_COUNT} V=1 VERBOSE=1 install -C bin

# install system configuration files
make -j ${CPU_COUNT} V=1 VERBOSE=1 install-sysconfDATA
