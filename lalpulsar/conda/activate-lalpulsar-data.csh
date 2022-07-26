#!/usr/bin/env csh
#
# Configure a conda environment for LALPulsar
#

# backup the environment's current setting
if ($?LALPULSAR_DATADIR) then
	setenv CONDA_BACKUP_LALPULSAR_DATADIR "${LALPULSAR_DATADIR}"
endif

# set the variable
setenv LALPULSAR_DATADIR "${CONDA_PREFIX}/share/lalpulsar"
