#!/usr/bin/env csh
#
# Deconfigure a conda environment for LALPulsar
#

if ($?LALPULSAR_DATADIR) then
	setenv LALPULSAR_DATADIR "$CONDA_BACKUP_LALPULSAR_DATADIR"
	unsetenv CONDA_BACKUP_LALPULSAR_DATADIR
else
	unsetenv LALPULSAR_DATADIR
endif
