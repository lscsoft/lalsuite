#!/usr/bin/env csh
#
# Deconfigure a conda environment for LALSimulation
#

if ($?LALSIMULATION_DATADIR) then
	setenv LALSIMULATION_DATADIR "$CONDA_BACKUP_LALSIMULATION_DATADIR"
	unsetenv CONDA_BACKUP_LALSIMULATION_DATADIR
else
	unsetenv LALSIMULATION_DATADIR
endif
