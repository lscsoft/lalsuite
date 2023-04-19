#!/usr/bin/env csh
#
# Configure a conda environment for LALSimulation
#

# backup the environment's current setting
if ($?LALSIMULATION_DATADIR) then
	setenv CONDA_BACKUP_LALSIMULATION_DATADIR "${LALSIMULATION_DATADIR}"
endif

# set the variable
setenv LALSIMULATION_DATADIR "${CONDA_PREFIX}/share/lalsimulation"
