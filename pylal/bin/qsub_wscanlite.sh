#!/bin/sh
if [ $# -ne 3 ]; then
  echo "Some of the arguments are missing, arguments needed are :"
  echo " - event GPS time"
  echo " - configuration file"
  echo " - output directory"
  exit
fi
eventTime=$1
curDirectory=`pwd`
if [ -e ./$2 ]; then
   configFile="$curDirectory/$2"
else
   configFile=$2
fi
#if [ -d ./$3 ]; then
#   outputDirectory="$curDirectory/$3"
#else
#   outputDirectory=$3
#fi
outputDirectory="$curDirectory/$3"

echo "Submitting : $curDirectory/SCRIPTS/wscanlite_in2p3.sh $eventTime $configFile $outputDirectory"

# This is how the job were set up before the switch to SL5.
#qsub -l platform=LINUX,u_sps_virgo,u_hpss_virgo,u_xrootd,matlab -q T -l M=2048MB <<eof 

qsub -l platform=LINUX,T=300000,scratch=1GB,M=3000MB,u_sps_virgo,u_hpss_virgo,matlab <<eof
$curDirectory/SCRIPTS/wscanlite_in2p3.sh $eventTime $configFile $outputDirectory
eof
