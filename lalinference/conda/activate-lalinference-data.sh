#!/usr/bin/env bash
#
# Configure a conda environment for LALInference
#

# preserve the user's existing setting
if [ ! -z "${LALINFERENCE_DATADIR+x}" ]; then
	export CONDA_BACKUP_LALINFERENCE_DATADIR="${LALINFERENCE_DATADIR}"
fi

# set the variable
export LALINFERENCE_DATADIR="${CONDA_PREFIX}/share/lalinference"
