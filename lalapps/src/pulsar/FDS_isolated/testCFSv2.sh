#!/bin/sh

## take user-arguments for CFS-v2:
extra_args="$@"

## allow 'make test' to work from builddir != srcdir
if [ -n "${srcdir}" ]; then
    builddir="./";
    injectdir="../Injections/"
else
    srcdir=.
fi

##---------- names of codes and input/output files
mfd_code="${injectdir}lalapps_Makefakedata_v4"
saf_code="${builddir}lalapps_SemiAnalyticF"
cfs_code="${builddir}lalapps_ComputeFStatistic"
cfsv2_code="${builddir}lalapps_ComputeFStatistic_v2"
cmp_code="${builddir}lalapps_compareFstats"

SFTdir="./testSFTs"

# test if LAL_DATA_PATH has been set ... needed to locate ephemeris-files
if [ -z "$LAL_DATA_PATH" ]; then
    if [ -n "$LAL_PREFIX" ]; then
	export LAL_DATA_PATH=".:${LAL_PREFIX}/share/lal";
    else
	echo
	echo "Need environment-variable LAL_PREFIX, or LAL_DATA_PATH to be set"
	echo "to your ephemeris-directory (e.g. /usr/local/share/lal)"
	echo "This might indicate an incomplete LAL installation"
	echo
	exit 1
    fi
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
mfd_CL="${saf_CL} --fmin=$mfd_fmin --Band=$mfd_FreqBand --Freq=$Freq --outSFTbname=$SFTdir/testSFT --f1dot=$f1dot --refTime=$refTime --outSFTv1"
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
echo "STEP 2: run CFS_v1 with perfect match"
echo "----------------------------------------------------------------------"
echo
outfile_v1="Fstat_v1.dat";
## common cmdline-options for v1 and v2
cfs_CL="--IFO=$IFO --Freq=$Freq --Alpha=$Alpha --Delta=$Delta --f1dot=$f1dot --DataFiles='$SFTdir/testSFT*' --refTime=$refTime"
if [ "$haveNoise" = false ]; then
    cfs_CL="$cfs_CL --SignalOnly"
fi

cmdline="$cfs_code $cfs_CL  --outputFstat=$outfile_v1 --expLALDemod=0 --Fthreshold=0";
echo $cmdline;

if ! eval $cmdline; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo " STEP 3: run CFS_v2 with perfect match"
echo "----------------------------------------------------------------------"
echo
outfile_v2NWon="Fstat_v2NWon.dat";
cmdlineNoiseWeightsOn="$cfsv2_code $cfs_CL --outputFstat=$outfile_v2NWon --TwoFthreshold=0 --UseNoiseWeights=true $extra_args";
echo $cmdlineNoiseWeightsOn;
if ! eval $cmdlineNoiseWeightsOn; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1;
fi

outfile_v2NWoff="Fstat_v2NWoff.dat";
cmdlineNoiseWeightsOff="$cfsv2_code $cfs_CL --outputFstat=$outfile_v2NWoff --TwoFthreshold=0 --UseNoiseWeights=false $extra_args";
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
echo "----- CFS v1: "
cat $outfile_v1 | sed -e"/^%.*/{d}"
echo "----- CFS v2 WITH noise-weights: "
cat $outfile_v2NWon | sed -e"/^%.*/{d}"
echo "----- CFS v2 WITHOUT noise-weights: "
cat $outfile_v2NWoff | sed -e"/^%.*/{d}"


echo
cmdline="$cmp_code -1 ./$outfile_v1 -2 ./$outfile_v2NWoff --clusterFiles=0 --Ftolerance=$Ftolerance"
echo $cmdline
if ! eval $cmdline; then
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi

cmdline="$cmp_code -1 ./$outfile_v1 -2 ./$outfile_v2NWon --clusterFiles=0 --Ftolerance=$Ftolerance"
echo $cmdline
if ! eval $cmdline; then
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi

## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm -rf $SFTdir $outfile_v1 $outfile_v2NWon $outfile_v2NWoff Fstats Fstats.log
fi
