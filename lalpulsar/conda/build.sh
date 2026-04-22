#!/bin/bash

set -ex

# load common options
. ${RECIPE_DIR}/common.sh

# use out-of-tree build
mkdir -pv _build
cd _build

# path to ephemeris files in `solar_system_ephemerides` package
SSE='./python*/site-packages/solar_system_ephemerides/ephemerides'

# configure
${SRC_DIR}/configure \
  ${CONFIGURE_ARGS} \
  --disable-python \
  --disable-swig-python \
  --enable-swig-iface \
  --with-fallback-data-path="\$(pkgdatadir):${SSE}/earth:${SSE}/sun" \
;

# build
${_make}

# test
if [[ "${CONDA_BUILD_CROSS_COMPILATION:-}" != "1"  ]]; then
  ${_make} check
fi
