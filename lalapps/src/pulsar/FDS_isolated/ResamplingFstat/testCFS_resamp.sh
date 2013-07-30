#!/bin/sh

## run all LALApps programs with memory debugging
export LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},memdbg"

## take user-arguments for CFS-resamp:
extra_args="$@"

builddir="./";
injectdir="../../Injections/"

##---------- names of codes and input/output files
mfd_code="${injectdir}lalapps_Makefakedata_v4"
saf_code="${builddir}../lalapps_SemiAnalyticF"
cfs2_code="${builddir}../lalapps_ComputeFStatistic_v2"
cfs_resamp_code="${builddir}lalapps_ComputeFStatistic_resamp"
cmp_code="${builddir}../lalapps_compareFstats"

SFTdir="./testCFSv2_resamp_sfts"

if [ -n "${LALPULSAR_DATADIR}" ]; then
    mfd_code="${mfd_code} -E ${LALPULSAR_DATADIR}"
    saf_code="${saf_code} -E ${LALPULSAR_DATADIR}"
    cfs_resamp_code="${cfs_resamp_code} -E ${LALPULSAR_DATADIR}"
fi

if [ -z "${LAL_DATA_PATH}" ]; then
    echo
    echo "Need environment-variable LAL_DATA_PATH to be set to include"
    echo "your ephemeris-directory (e.g. /usr/local/share/lalpulsar)"
    echo "This might indicate an incomplete LAL+LALPULSAR installation"
    echo
    exit 1
fi

Ftolerance=0.05
# ---------- fixed parameter of our test-signal
Tsft=1800;
startTime=711595934
refTime=701595833  ## $startTime
duration=144000		## 40 hours

mfd_FreqBand=2.0;

Alpha=2.0
Delta=-0.5

h0=1
cosi=-0.3

psi=0.6
phi0=1.5

Freq=100.12345
cfsFreqBand=4e-6;
dFreq=6.944444444e-06

mfd_fmin=$(echo $Freq $mfd_FreqBand | awk '{printf "%g", $1 - $2 / 2.0}');

f1dot=-1e-10;

noiseSqrtSh=5

## ------------------------------------------------------------

if [ "$noiseSqrtSh" != 0 ]; then
    sqrtSh=$noiseSqrtSh
    haveNoise=true;
else
    sqrtSh=1;	## for SemiAnalyticF signal-only case
    haveNoise=false;
fi

IFO=LHO

##--------------------------------------------------
## test starts here
##--------------------------------------------------

echo
echo "----------------------------------------------------------------------"
echo " STEP 1: Generate Fake Signal"
echo "----------------------------------------------------------------------"
echo
if [ ! -d "$SFTdir" ]; then
    mkdir $SFTdir;
else
    rm -f $SFTdir/*;
fi

# this part of the command-line is compatible with SemiAnalyticF:
saf_CL=" --Alpha=$Alpha --Delta=$Delta --IFO=$IFO --Tsft=$Tsft --startTime=$startTime --duration=$duration --h0=$h0 --cosi=$cosi --psi=$psi --phi0=$phi0"
# concatenate this with the mfd-specific switches:
mfd_CL="${saf_CL} --fmin=$mfd_fmin --Band=$mfd_FreqBand --Freq=$Freq --outSFTbname=$SFTdir --f1dot=$f1dot --refTime=$refTime"
if [ "$haveNoise" = true ]; then
    mfd_CL="$mfd_CL --noiseSqrtSh=$sqrtSh";
fi

cmdline="$mfd_code $mfd_CL";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

echo
echo -n "Running '$saf_code' ... "
cmdline="$saf_code $saf_CL --sqrtSh=$sqrtSh"
echo $cmdline
if ! resF=`eval $cmdline 2> /dev/null`; then
    echo "Error ... something failed running '$saf_code' ..."
    exit 1;
fi
echo  "ok."
res2F=`echo $resF | awk '{printf "%g", 2.0 * $1}'`
echo "The SemiAnalyticF calculations predicts: 2F = $res2F"

echo
echo "----------------------------------------------------------------------"
echo "STEP 2: run CFS_v2 with perfect match"
echo "----------------------------------------------------------------------"
echo
outfile_v2="Fstat_v2.dat";
## common cmdline-options for v1 and resamp
cfs_CL="--IFO=$IFO --Freq=$Freq --FreqBand=$cfsFreqBand --dFreq=$dFreq --Alpha=$Alpha --Delta=$Delta --f1dot=$f1dot --DataFiles='${SFTdir}/*.sft' --refTime=$refTime"
if [ "$haveNoise" = false ]; then
    cfs_CL="$cfs_CL --SignalOnly"
fi

cmdline="$cfs2_code $cfs_CL  --outputFstat=$outfile_v2 --TwoFthreshold=0";
echo $cmdline;

if ! eval $cmdline; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo " STEP 3: run CFS_resamp with perfect match"
echo "----------------------------------------------------------------------"
echo
outfile_resampNWon="Fstat_resampNWon.dat";
cmdlineNoiseWeightsOn="$cfs_resamp_code $cfs_CL --outputFstat=$outfile_resampNWon --TwoFthreshold=0 --UseNoiseWeights=true $extra_args";
echo $cmdlineNoiseWeightsOn;
if ! eval $cmdlineNoiseWeightsOn; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1;
fi

outfile_resampNWoff="Fstat_resampNWoff.dat";
cmdlineNoiseWeightsOff="$cfs2_code $cfs_CL --outputFstat=$outfile_resampNWoff --TwoFthreshold=0 --UseNoiseWeights=false $extra_args";
echo $cmdlineNoiseWeightsOff;
if ! eval $cmdlineNoiseWeightsOff; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1;
fi

echo
echo "----------------------------------------"
echo " STEP 4: Comparing results: "
echo "----------------------------------------"
echo
echo "----- CFS v2: "
cat $outfile_v2 | grep -v "%%"
echo "----- CFS resamp WITH noise-weights: "
cat $outfile_resampNWon | grep -v "%%"
echo "----- CFS resamp WITHOUT noise-weights: "
cat $outfile_resampNWoff | grep -v "%%"


echo
cmdline="$cmp_code -1 ./$outfile_v2 -2 ./$outfile_resampNWoff --clusterFiles=0 --Ftolerance=$Ftolerance"
echo $cmdline
if ! eval $cmdline; then
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi

cmdline="$cmp_code -1 ./$outfile_v2 -2 ./$outfile_resampNWon --clusterFiles=0 --Ftolerance=$Ftolerance"
echo $cmdline
if ! eval $cmdline; then
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi

## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm -rf $SFTdir $outfile_v2 $outfile_resampNWon $outfile_resampNWoff Fstats Fstats.log
fi
