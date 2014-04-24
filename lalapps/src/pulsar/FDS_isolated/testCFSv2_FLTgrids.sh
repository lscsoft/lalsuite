#!/bin/sh

## set LAL debug level
echo "Setting LAL_DEBUG_LEVEL=${LAL_DEBUG_LEVEL:-msglvl1,memdbg}"
export LAL_DEBUG_LEVEL

## allow 'make test' to work from builddir != srcdir
if [ -z "${srcdir}" ]; then
    srcdir=`dirname $0`
fi

builddir="./";
injectdir="../Injections/"

## user-controlled level of debug-output detail
if [ -n "$DEBUG" ]; then
    debug=${DEBUG}
else
    debug=0	## default=quiet
fi

## names of codes and input/output files
mfd_code="${injectdir}lalapps_Makefakedata_v5"
cfsv2_code="${builddir}lalapps_ComputeFStatistic_v2"
cmp_code="${builddir}lalapps_compareFstats"

## check for ephemeris path
if [ -z "${LAL_DATA_PATH}" ]; then
    echo
    echo "Need environment-variable LAL_DATA_PATH to be set to include"
    echo "your ephemeris-directory (e.g. /usr/local/share/lalpulsar)"
    echo "This might indicate an incomplete LAL+LALPULSAR installation"
    echo
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo " STEP 1: Generate Fake SFTs"
echo "----------------------------------------------------------------------"
echo

## create SFT directory
SFTdir="./testCFSv2_FLTgrids"
if [ ! -d "$SFTdir" ]; then
    mkdir $SFTdir;
else
    rm -f $SFTdir/*;
fi

## create timestamps file
cat <<EOF >"${SFTdir}/timestamps.txt"
800000000 0
800432000 0
EOF

## build MFD_v5 command line
mfd_CL="--outSingleSFT=true --outSFTdir=${SFTdir} --IFOs=H1 --sqrtSX=1.0 --timestampsFiles=${SFTdir}/timestamps.txt --fmin=100 --Band=1.0 --randSeed=12345"

## generate SFTs
cmdline="$mfd_code $mfd_CL";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error: something failed when running '$mfd_code' ..."
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo " STEP 2: Run CFS_v2 with various FlatLatticeTiling-generated grid types"
echo "----------------------------------------------------------------------"
echo

## build CFS_v2 command line for --gridType=8
gridFile8="./testCFSv2_FLTgrids_8.dat"
cfsv2_CL="--Alpha=6.1 --Delta=1.2 --Freq=100.4 --FreqBand=5e-4 --f1dot=-1e-10 --f1dotBand=1e-10 --DataFiles='${SFTdir}/*.sft' --TwoFthreshold=0 --gridType=8 --metricMismatch=0.5 --outputFstat=${gridFile8}"

## run CFS_v2
cmdline="$cfsv2_code $cfsv2_CL";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error: something failed when running '$cfsv2_code' ..."
    exit 1
fi

## build CFS_v2 command line for --gridType=9
gridFile9="./testCFSv2_FLTgrids_9.dat"
cfsv2_CL="--Alpha=6.1 --Delta=1.2 --Freq=100.4 --FreqBand=8e-5 --spindownAge=1e11 --minBraking=2 --maxBraking=5 --DataFiles='${SFTdir}/*.sft' --TwoFthreshold=0 --gridType=9 --metricMismatch=0.5 --outputFstat=${gridFile9}"

## run CFS_v2
cmdline="$cfsv2_code $cfsv2_CL";
echo $cmdline;
if ! eval $cmdline; then
   echo "Error: something failed when running '$cfsv2_code' ..."
   exit 1
fi

echo
echo "----------------------------------------"
echo " STEP 3: Comparing results"
echo "----------------------------------------"
echo

## compare --gridType=8
echo "Comparing gridType=8:"
cmdline="$cmp_code -1 ${gridFile8} -2 ${srcdir}/${gridFile8}.ref.gz --clusterFiles=0 --Ftolerance=0.1";
echo $cmdline
if ! eval $cmdline; then
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi

## compare --gridType=9
echo "Comparing gridType=9:"
cmdline="$cmp_code -1 ${gridFile9} -2 ${srcdir}/${gridFile9}.ref.gz --clusterFiles=0 --Ftolerance=0.1";
echo $cmdline
if ! eval $cmdline; then
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi

echo
echo "----------------------------------------"
echo " DONE!"
echo "----------------------------------------"
echo

## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm -rf ${SFTdir} ${gridFile8} ${gridFile9}
fi
