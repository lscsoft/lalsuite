#!/usr/bin/env csh
#
# Configure a conda environment for LALInference
#

# backup the environment's current setting
if ($?LALINFERENCE_DATADIR) then
	setenv CONDA_BACKUP_LALINFERENCE_DATADIR "${LALINFERENCE_DATADIR}"
endif

# set the variable
setenv LALINFERENCE_DATADIR "${CONDA_PREFIX}/share/lalinference"
