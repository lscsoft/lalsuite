#!/bin/bash 

# source the .profile 
source ${HOME}/.profile 
echo '... sourced .profile on node' $1 

# launch the python script 
echo '... launching the job on node ' $1
python /home/anuradha.gupta/BoxingDay/Prec_Fit/condor_runs/GW151226/run_condor_output_evolved_spins.py GW151226_IMRPp_SpinFix /home/anuradha.gupta/BoxingDay/Prec_Fit/condor_runs/GW151226 $1
