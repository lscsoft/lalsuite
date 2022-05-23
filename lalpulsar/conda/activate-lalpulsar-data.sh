#!/usr/bin/env bash
#
# Configure a conda environment for LALPulsar
#

# preserve the user's existing setting
if [ ! -z "${LALPULSAR_DATADIR+x}" ]; then
	export CONDA_BACKUP_LALPULSAR_DATADIR="${LALPULSAR_DATADIR}"
fi

# set the variable
export LALPULSAR_DATADIR="${CONDA_PREFIX}/share/lalpulsar"
