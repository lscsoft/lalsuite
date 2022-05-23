#!/usr/bin/env bash
#
# De-configure a conda environment for LALPulsar
#

# restore from backup
if [ ! -z "${CONDA_BACKUP_LALPULSAR_DATADIR}" ]; then
	export LALPULSAR_DATADIR="${CONDA_BACKUP_LALPULSAR_DATADIR}"
	unset CONDA_BACKUP_LALPULSAR_DATADIR
# no backup, just unset
else
	unset LALPULSAR_DATADIR
fi
