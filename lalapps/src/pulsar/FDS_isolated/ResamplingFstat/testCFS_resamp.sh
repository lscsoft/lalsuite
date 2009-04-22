#!/bin/sh

## take user-arguments for CFS
extra_args="$@"

##---------- names of codes and input/output files
saf_code="lalapps_SemiAnalyticF"
mfd_code="lalapps_Makefakedata_v4"
cfs_code="lalapps_ComputeFStatistic"
cfs_resamp_code="lalapps_ComputeFStatistic_resamp"
cmp_code="lalapps_compareFstats"

SFTdir="./_testSFTs"

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
cfsFreqBand=4e-5;
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
outfile_v1="./Fstat_v1.dat";
## common cmdline-options for v1 and _resamp
cfs_CL="--IFO=$IFO --Freq=$Freq --FreqBand=$cfsFreqBand --Alpha=$Alpha --Delta=$Delta --f1dot=$f1dot --DataFiles='$SFTdir/testSFT*' --refTime=$refTime"
if [ "$haveNoise" = false ]; then
    cfs_CL="$cfs_CL --SignalOnly"
fi

cmdline="$cfs_code $cfs_CL  --outputFstat=$outfile_v1 --expLALDemod=0 --Fthreshold=0 --dFreq=$dFreq";
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
outfile_resamp="./Fstat_resamp.dat";
cmdline_resamp="$cfs_resamp_code $cfs_CL --outputFstat=$outfile_resamp --TwoFthreshold=0 --UseNoiseWeights=false $extra_args";
echo $cmdline_resamp;
if ! eval $cmdline_resamp; then
    echo "Error.. something failed when running '$cfs_resamp_code' ..."
    exit 1;
fi

echo
echo "----------------------------------------"
echo " STEP 4: Comparing results: "
echo "----------------------------------------"
echo
echo "----- CFS v1: "
cat $outfile_v1 | sed -e"/^%.*/{d}"
echo "----- CFS resamp (with noise-weights turned OFF): "
cat $outfile_resamp | sed -e"/^%.*/{d}"


echo
cmdline="$cmp_code -1 $outfile_v1 -2 $outfile_resamp --clusterFiles=0 --Ftolerance=$Ftolerance"
echo $cmdline
if ! eval $cmdline; then
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi
