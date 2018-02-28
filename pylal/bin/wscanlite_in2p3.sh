#!/bin/sh
#PBS -V
#PBS -p u=99
#PBS -l platform=LINUX,T=300000,scratch=1GB,M=3000MB,u_sps_virgo,u_hpss_virgo,matlab

# CCIN2P3 specific script calling omega scan to analyze events in VSR1
# Call with
#      wscan_ccin2p3.sh event_time config_file output_directory
#
# where event time is the GPS time of the event in seconds
#       config_file is the qscan configuration file
#       output directory is the directory for outputing the html result
#
#     example : qscan-ccin2p3.sh 863559000 V1-raw-gravitational.txt test
#

if [ $# -lt 3 ]; then
  echo "Some of the arguments are missing, arguments needed are :"
  echo " - event GPS time"
  echo " - configuration file"
  echo " - output directory"
  exit
fi

eventTime=$1
configFile=$2
outputDirectory=$3

omegaDirectory="$THRONG_DIR/pro/omegadev/omega_r2757"
matlabDirectory="/afs/in2p3.fr/system/amd64_sl5/usr/local/matlabR2009a/bin/glnxa64"
FFLFile="/afs/in2p3.fr/group/virgo/BKDB/VSR2/VSR2_raw.ffl"

# Set path for omega
testpath=`echo $PATH | grep -i 'omegadev/omega_r2757/bin'`

if [ -z $testpath ]; then
  export PATH=$omegaDirectory/bin:$PATH
fi

# Set ld_library_path for matlab used by omega
testpath=`echo $PATH | grep -i 'matlabR2009a/bin/glnxa64'`

if [ -z $testpath ]; then
  export LD_LIBRARY_PATH=$matlabDirectory:$LD_LIBRARY_PATH
fi

# LOAD lscsoft
source /afs/in2p3.fr/throng/virgo/pro/lscsoft/lscsoft-user-env.sh

if [ -d $outputDirectory ]; then
  echo ""
  echo "Directory $outputDirectory exists already."
  echo "**********************"
  echo "*** Cleaning it... ***"
  echo "**********************"
  rm -rf $outputDirectory
fi

# Execute the wscanlite

OMEGASCAN="$omegaDirectory/bin/wpipeline scan -c $configFile -f $FFLFile -o $outputDirectory $eventTime"

echo "execute : $OMEGASCAN"
export LD_LIBRARY_PATH_SAV=${LD_LIBRARY_PATH}
source /usr/local/shared/bin/xrootd_env.sh
$OMEGASCAN
unset LD_PRELOAD
unset XROOTD_VMP
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH_SAV}

exit 0
