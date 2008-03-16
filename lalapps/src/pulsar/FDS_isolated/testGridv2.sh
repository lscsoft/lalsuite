#!/bin/sh

## take user-arguments for CFS-v2:
extra_args="$@"
debug=0;

## allow 'make test' to work from builddir != srcdir
if [ -n "${srcdir}" ]; then
    builddir="./";
    injectdir="../Injections/"
else
    srcdir=.
fi

##---------- names of codes and input/output files
mfd_code="${injectdir}lalapps_Makefakedata_v4"
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

# ---------- fixed parameter of our test-signal
Tsft=1800;
startTime=711595934
refTime=701595833  ## $startTime
duration=144000		## 40 hours

Alpha=2.0
AlphaBand=0.1
dAlpha=0.04
Delta=-0.5
DeltaBand=0.1
dDelta=0.04

h0=1
cosi=-0.3

psi=0.6
phi0=1.5

Freq=100.12345
FreqBand=0.1
dFreq=0.04
f1dot=-1e-13
f1dotBand=1e-14
df1dot=4e-15

mfd_FreqBand=4.0;
mfd_fmin=$(echo $Freq $mfd_FreqBand | awk '{printf "%g", $1 - $2 / 2.0}');

echo "mfd_fmin = $mfd_fmin"

noiseSqrtSh=0

## ------------------------------------------------------------

if [ "$noiseSqrtSh" != 0 ]; then
    sqrtSh=$noiseSqrtSh
    haveNoise=true;
else
    sqrtSh=1;	## for SemiAnalyticF signal-only case
    haveNoise=false;
fi

IFO=LLO

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
##echo -n "Running '$saf_code' ... "
##cmdline="$saf_code $saf_CL --sqrtSh=$sqrtSh"
##echo $cmdline
##if ! resF=`eval $cmdline 2> /dev/null`; then
##    echo "Error ... something failed running '$saf_code' ..."
##    exit 1;
##fi
##echo  "ok."
##res2F=`echo $resF | awk '{printf "%g", 2.0 * $1}'`
##echo "The SemiAnalyticF calculations predicts: 2F = $res2F"

echo
echo "----------------------------------------------------------------------"
echo "STEP 2: run CFS_v1 with various grid-type"
echo "----------------------------------------------------------------------"
echo
gridType0=" --gridType=0";	## flat grid
gridType1=" --gridType=1";	## isotropic
gridType2=" --gridType=2 --metricType=1 --metricMismatch=0.1";	## metric grid
outputv1_0="./Fstatv1_grid0.dat";
outputv1_1="./Fstatv1_grid1.dat";
outputv1_2="./Fstatv1_grid2.dat";

## common cmdline-options for v1 and v2
cfs_CL="--IFO=$IFO --Freq=$Freq --FreqBand=$FreqBand --dFreq=$dFreq --Alpha=$Alpha --AlphaBand=$AlphaBand --dAlpha=$dAlpha --Delta=$Delta --DeltaBand=$DeltaBand --dDelta=$dDelta --f1dot=$f1dot --f1dotBand=$f1dotBand --df1dot=$df1dot --DataFiles='$SFTdir/testSFT*' --refTime=$refTime"
if [ "$haveNoise" = false ]; then
    cfs_CL="$cfs_CL --SignalOnly"
fi

cmdlinev1="$cfs_code $cfs_CL --expLALDemod=1 --Fthreshold=0"

## ----- grid=0
echo -n "CFSv1 using gridType=0 ..."
cmd0="$cmdlinev1 $gridType0 --outputFstat=$outputv1_0";
if [ $debug = 1 ]; then echo; echo $cmd0; fi
if ! eval $cmd0; then
    echo "Error.. something failed when running '$cmd0' ..."
    exit 1
fi
echo "done."

## ----- grid=1
echo -n "CFSv1 using gridType=1 ..."
cmd0="$cmdlinev1 $gridType1 --outputFstat=$outputv1_1";
if [ $debug = 1 ]; then echo; echo $cmd0; fi
if ! eval $cmd0; then
    echo "Error.. something failed when running '$cmd1' ..."
    exit 1
fi
echo "done."

## ----- grid=2
echo -n "CFSv1 using gridType=2 ..."
cmd0="$cmdlinev1 $gridType2 --outputFstat=$outputv1_2";
if [ $debug = 1 ]; then echo; echo $cmd0; fi
if ! eval $cmd0; then
    echo "Error.. something failed when running '$cmd2' ..."
    exit 1
fi
echo "done."


echo
echo "----------------------------------------------------------------------"
echo " STEP 3: run CFS_v2 for various grid-types"
echo "----------------------------------------------------------------------"
echo
outputv2_0="./Fstatv2_grid0.dat";
outputv2_1="./Fstatv2_grid1.dat";
outputv2_2="./Fstatv2_grid2.dat";
cmdlinev2="$cfsv2_code $cfs_CL --TwoFthreshold=0 $extra_args";

## ----- grid=0
echo -n "CFSv2 using gridType=0 ..."
cmd0="$cmdlinev2 $gridType0 --outputFstat=$outputv2_0";
if [ $debug = 1 ]; then echo; echo $cmd0; fi
if ! eval $cmd0; then
    echo "Error.. something failed when running '$cmd0' ..."
    exit 1
fi
echo "done."

## ----- grid=1
echo -n "CFSv2 using gridType=1 ..."
cmd0="$cmdlinev2 $gridType1 --outputFstat=$outputv2_1";
if [ $debug = 1 ]; then echo; echo $cmd0; fi
if ! eval $cmd0; then
    echo "Error.. something failed when running '$cmd1' ..."
    exit 1
fi
echo "done."

## ----- grid=2
echo -n "CFSv2 using gridType=2 ..."
cmd0="$cmdlinev2 $gridType2 --outputFstat=$outputv2_2";
if [ $debug = 1 ]; then echo; echo $cmd0; fi
if ! eval $cmd0; then
    echo "Error.. something failed when running '$cmd2' ..."
    exit 1
fi
echo "done."


echo
echo "----------------------------------------"
echo " STEP 4: Comparing results: "
echo "----------------------------------------"
echo

## ----- grid=0
echo -n "Comparing gridType=0 ... "
if [ $debug = 1 ]; then echo; echo $cmd0; fi
cmd0="$cmp_code -1 ./$outputv1_0  -2 ./$outputv2_0 --clusterFiles=0 --Ftolerance=0.1";
if ! eval $cmd0; then
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi

## ----- grid=1
echo -n "Comparing gridType=1 ... "
if [ $debug = 1 ]; then echo; echo $cmd0; fi
cmd0="$cmp_code -1 ./$outputv1_1  -2 ./$outputv2_1 --clusterFiles=0 --Ftolerance=0.1";
if ! eval $cmd0; then
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi


## ----- grid=2
echo -n "Comparing gridType=2 ... "
if [ $debug = 1 ]; then echo; echo $cmd0; fi
cmd0="$cmp_code -1 ./$outputv1_2  -2 ./$outputv2_2 --clusterFiles=0 --Ftolerance=0.1";
if ! eval $cmd0; then
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi

echo

## clean up files
rm -rf $SFTdir $outputv2_0 $outputv2_1 $outputv2_2 $outputv1_0 $outputv1_1 $outputv1_2 Fstats Fstats.log

