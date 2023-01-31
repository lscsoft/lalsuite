#!/bin/bash
#
# Configure, built, test, and install the Python language bindings
# for a LALSuite subpackage.
#

set -e

# load common options
. ${RECIPE_DIR}/common.sh

# build python in a sub-directory using a copy of the C build
_builddir="_build${PY_VER}"
cp -r _build ${_builddir}
cd ${_builddir}

# configure only python bindings and pure-python extras
${SRC_DIR}/configure \
  ${CONFIGURE_ARGS} \
  --disable-swig-iface \
  --enable-sistr \
  --enable-python \
  --enable-swig-python \
;

# patch out dependency_libs from libtool archive to prevent overlinking
sed -i.tmp '/^dependency_libs/d' lib/lib${PKG_NAME##*-}.la

# build
${_make} -C swig LIBS=""
${_make} -C python LIBS=""

# test
if [[ "${CONDA_BUILD_CROSS_COMPILATION:-}" != "1" || "${CROSSCOMPILING_EMULATOR}" != "" ]]; then
  ${_make} check -C swig
fi

# install
${_make} -C swig install-exec  # swig bindings
${_make} -C python install  # pure-python extras
