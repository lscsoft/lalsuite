#!/bin/bash

set -ex

# load common options
. ${RECIPE_DIR}/common.sh

# build in a sub-directory using a copy of the python build
_builddir="_build${PY_VER}_${mpi}"
cp -r _build${PY_VER} ${_builddir}
cd ${_builddir}

# enable cross compiling with openmpi
if [[ "${CONDA_BUILD_CROSS_COMPILATION:-}" == "1" && "${mpi}" == "openmpi" ]]; then
  cp -rf ${PREFIX}/share/openmpi/*.txt ${BUILD_PREFIX}/share/openmpi/
fi

if [[ "${mpi}" != "nompi" ]]; then
  CONFIGURE_ARGS="${CONFIGURE_ARGS} --enable-mpi"
fi

# configure only what we need to build executables (and man pages)
${SRC_DIR}/configure \
  ${CONFIGURE_ARGS} \
  --disable-swig \
  --enable-python \
;

# install binaries
${_make} install -C bin

# install system configuration files
${_make} install-sysconfDATA
