#!/usr/bin/env bash
#
# De-configure a conda environment for LALInference
#

# restore from backup
if [ ! -z "${CONDA_BACKUP_LALINFERENCE_DATADIR}" ]; then
	export LALINFERENCE_DATADIR="${CONDA_BACKUP_LALINFERENCE_DATADIR}"
	unset CONDA_BACKUP_LALINFERENCE_DATADIR
# no backup, just unset
else
	unset LALINFERENCE_DATADIR
fi
