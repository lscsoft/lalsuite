#!/usr/bin/env csh
#
# Deconfigure a conda environment for LALInference
#

if ($?LALINFERENCE_DATADIR) then
	setenv LALINFERENCE_DATADIR "$CONDA_BACKUP_LALINFERENCE_DATADIR"
	unsetenv CONDA_BACKUP_LALINFERENCE_DATADIR
else
	unsetenv LALINFERENCE_DATADIR
endif
