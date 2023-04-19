#!/bin/bash

set -ex

# load common options
. ${RECIPE_DIR}/common.sh

# use out-of-tree build
mkdir -pv _build
cd _build

# configure
${SRC_DIR}/configure \
  ${CONFIGURE_ARGS} \
  --enable-python \
  --enable-swig-iface \
  --enable-swig-python \
;

# build the C library
${_make} -C include
${_make} -C lib

# patch out dependency_libs from libtool archive to prevent overlinking
sed -i.tmp '/^dependency_libs/d' lib/liblalsimulation.la

# build the SWIG bindings and Python modules
${_make} -C swig LIBS=""
${_make} -C python LIBS=""

# build everything else
${_make}

# test
if [[ "${CONDA_BUILD_CROSS_COMPILATION:-}" != "1" || "${CROSSCOMPILING_EMULATOR}" != "" ]]; then
  ${_make} check
fi

# install
${_make} install

# install activate/deactivate scripts
for action in activate deactivate; do
	mkdir -p ${PREFIX}/etc/conda/${action}.d
	for ext in sh csh; do
		_target="${PREFIX}/etc/conda/${action}.d/${action}-lalsimulation.${ext}"
		echo "-- Installing: ${_target}"
		cp "${RECIPE_DIR}/${action}-lalsimulation.${ext}" "${_target}"
	done
done
