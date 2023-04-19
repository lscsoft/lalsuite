#!/usr/bin/env bash
#
# De-configure a conda environment for LALSimulation
#

# restore from backup
if [ ! -z "${CONDA_BACKUP_LALSIMULATION_DATADIR}" ]; then
	export LALSIMULATION_DATADIR="${CONDA_BACKUP_LALSIMULATION_DATADIR}"
	unset CONDA_BACKUP_LALSIMULATION_DATADIR
# no backup, just unset
else
	unset LALSIMULATION_DATADIR
fi
