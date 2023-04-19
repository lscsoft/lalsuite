#!/usr/bin/env bash
#
# Configure a conda environment for LALSimulation
#

# preserve the user's existing setting
if [ ! -z "${LALSIMULATION_DATADIR+x}" ]; then
	export CONDA_BACKUP_LALSIMULATION_DATADIR="${LALSIMULATION_DATADIR}"
fi

# set the variable
export LALSIMULATION_DATADIR="${CONDA_PREFIX}/share/lalsimulation"
