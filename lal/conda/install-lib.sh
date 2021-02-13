#!/bin/bash

set -ex

cd _build

# install library and headers
make -j ${CPU_COUNT} V=1 VERBOSE=1 -C lib install HDF5_LIBS="-L${PREFIX} -lhdf5 -lhdf5_hl"

# install SWIG binding definitions and headers
make -j ${CPU_COUNT} V=1 VERBOSE=1 -C swig install-data

# install pkg-config
make -j ${CPU_COUNT} V=1 VERBOSE=1 install-pkgconfigDATA

# -- create activate/deactivate scripts

# strip out the 'lib' package name prefix
LALSUITE_NAME=${PKG_NAME#"lib"}
LALSUITE_NAME_UPPER=$(echo ${LALSUITE_NAME} | awk '{ print toupper($0) }')

# activate.sh
ACTIVATE_SH="${PREFIX}/etc/conda/activate.d/activate_${PKG_NAME}.sh"
mkdir -p $(dirname ${ACTIVATE_SH})
cat > ${ACTIVATE_SH} << EOF
#!/bin/bash
export CONDA_BACKUP_${LALSUITE_NAME_UPPER}_DATADIR="\${${LALSUITE_NAME_UPPER}_DATADIR:-empty}"
export ${LALSUITE_NAME_UPPER}_DATADIR="/opt/anaconda1anaconda2anaconda3/share/${LALSUITE_NAME}"
EOF
# deactivate.sh
DEACTIVATE_SH="${PREFIX}/etc/conda/deactivate.d/deactivate_${PKG_NAME}.sh"
mkdir -p $(dirname ${DEACTIVATE_SH})
cat > ${DEACTIVATE_SH} << EOF
#!/bin/bash
if [ "\${CONDA_BACKUP_${LALSUITE_NAME_UPPER}_DATADIR}" = "empty" ]; then
	unset ${LALSUITE_NAME_UPPER}_DATADIR
else
	export ${LALSUITE_NAME_UPPER}_DATADIR="\${CONDA_BACKUP_${LALSUITE_NAME_UPPER}_DATADIR}"
fi
unset CONDA_BACKUP_${LALSUITE_NAME_UPPER}_DATADIR
EOF
