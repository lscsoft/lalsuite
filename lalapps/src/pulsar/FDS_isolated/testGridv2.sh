#!/bin/sh

## run all LALApps programs with memory debugging
export LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},memdbg"

## take user-arguments for CFS-v2:
extra_args="$@"

builddir="./";
injectdir="../Injections/"

## ----- user-controlled level of debug-output detail
if [ -n "$DEBUG" ]; then
    debug=${DEBUG}
else
    debug=0	## default=quiet
fi

##---------- names of codes and input/output files
mfd_code="${injectdir}lalapps_Makefakedata_v4"
cfs_code="${builddir}lalapps_ComputeFStatistic"
cfsv2_code="${builddir}lalapps_ComputeFStatistic_v2"
cmp_code="${builddir}lalapps_compareFstats"

SFTdir="./testGridv2_sfts"

if [ -n "${LALPULSAR_DATADIR}" ]; then
    mfd_code="${mfd_code} -E ${LALPULSAR_DATADIR}"
fi

if [ -z "${LAL_DATA_PATH}" ]; then
    echo
    echo "Need environment-variable LAL_DATA_PATH to be set to include"
    echo "your ephemeris-directory (e.g. /usr/local/share/lalpulsar)"
    echo "This might indicate an incomplete LAL+LALPULSAR installation"
    echo
    exit 1
fi

# ---------- fixed parameter of our test-signal
Tsft=1800;
startTime=711595934
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
mfd_CL="${saf_CL} --fmin=$mfd_fmin --Band=$mfd_FreqBand --Freq=$Freq --outSFTbname=$SFTdir/testSFT --f1dot=$f1dot --outSFTv1"
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
sky_CL="--Alpha=$Alpha --AlphaBand=$AlphaBand --dAlpha=$dAlpha --Delta=$Delta --DeltaBand=$DeltaBand --dDelta=$dDelta"
spin_CL="--Freq=$Freq --FreqBand=$FreqBand --dFreq=$dFreq --f1dot=$f1dot --f1dotBand=$f1dotBand --df1dot=$df1dot"
cfs_CL="--IFO=$IFO --DataFiles='$SFTdir/testSFT*'"
if [ "$haveNoise" = false ]; then
    cfs_CL="$cfs_CL --SignalOnly"
fi

cmdlinev1="$cfs_code $cfs_CL --expLALDemod=1 --Fthreshold=0 ${sky_CL} ${spin_CL}"

## ----- grid=0
echo "CFSv1 using gridType=0:"
cmd0="$cmdlinev1 $gridType0 --outputFstat=$outputv1_0";
echo $cmd0
if ! eval $cmd0; then
    echo "Error.. something failed when running '$cmd0' ..."
    exit 1
fi

## ----- grid=1
echo "CFSv1 using gridType=1:"
cmd0="$cmdlinev1 $gridType1 --outputFstat=$outputv1_1";
echo $cmd0
if ! eval $cmd0; then
    echo "Error.. something failed when running '$cmd0' ..."
    exit 1
fi

## ----- grid=2
echo "CFSv1 using gridType=2:"
cmd0="$cmdlinev1 $gridType2 --outputFstat=$outputv1_2";
echo $cmd0
if ! eval $cmd0; then
    echo "Error.. something failed when running '$cmd0' ..."
    exit 1
fi


echo
echo "----------------------------------------------------------------------"
echo " STEP 3: run CFS_v2 for various grid-types"
echo "----------------------------------------------------------------------"
echo
outputv2_0="./Fstatv2_grid0.dat";
outputv2_1="./Fstatv2_grid1.dat";
outputv2_2="./Fstatv2_grid2.dat";
outputv2_6="./Fstatv2_grid6.dat";

cmdlinev2="$cfsv2_code $cfs_CL --TwoFthreshold=0 --Dterms=16 $extra_args";

## ----- grid=0
echo "CFSv2 using gridType=0:"
cmd0="$cmdlinev2 ${sky_CL} ${spin_CL} $gridType0 --outputFstat=$outputv2_0";
echo $cmd0
if ! eval $cmd0; then
    echo "Error.. something failed when running '$cmd0' ..."
    exit 1
fi

## ----- grid=1
echo "CFSv2 using gridType=1:"
cmd0="$cmdlinev2 ${sky_CL} ${spin_CL} $gridType1 --outputFstat=$outputv2_1";
echo $cmd0
if ! eval $cmd0; then
    echo "Error.. something failed when running '$cmd0' ..."
    exit 1
fi

## ----- grid=2
echo "CFSv2 using gridType=2:"
cmd0="$cmdlinev2 ${sky_CL} ${spin_CL} $gridType2 --outputFstat=$outputv2_2";
echo $cmd0
if ! eval $cmd0; then
    echo "Error.. something failed when running '$cmd2' ..."
    exit 1
fi

## ----- recompute with --gridFile option --gridType=6
echo "recompute with CFSv2 using gridType=6:"
## extract a 6-column gridFile from the previous result output-file
gridFile="gridv2_2.dat";
rm -f ${gridFile};
awk_extract6='{printf "%s %s %s %s %s %s\n", $1, $2, $3, $4, $5, $6 >> "gridv2_2.dat" }'
grid_line=$(sed '/^%.*/d' ${outputv2_2} | awk "$awk_extract6")


cmd0="$cmdlinev2 --gridType=6 --gridFile=./${gridFile} --outputFstat=$outputv2_6";
echo $cmd0
if ! eval $cmd0; then
    echo "Error.. something failed when running '$cmd2' ..."
    exit 1
fi


echo
echo "----------------------------------------"
echo " STEP 4: Comparing results: "
echo "----------------------------------------"
echo

## ----- grid=0
echo "Comparing gridType=0:"
echo $cmd0
cmd0="$cmp_code -1 ./$outputv1_0  -2 ./$outputv2_0 --clusterFiles=0 --Ftolerance=0.1";
if ! eval $cmd0; then
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi

## ----- grid=1
echo "Comparing gridType=1:"
echo $cmd0
cmd0="$cmp_code -1 ./$outputv1_1  -2 ./$outputv2_1 --clusterFiles=0 --Ftolerance=0.1";
if ! eval $cmd0; then
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi


## ----- grid=2
echo "Comparing gridType=2:"
echo $cmd0
cmd0="$cmp_code -1 ./$outputv1_2  -2 ./$outputv2_2 --clusterFiles=0 --Ftolerance=0.1";
if ! eval $cmd0; then
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi

## ----- grid=6
echo "Comparing gridType=6:"
echo $cmd0
cmd0="$cmp_code -1 ./$outputv2_2 -2 ./$outputv2_6 --clusterFiles=0 --Ftolerance=0.01";
if ! eval $cmd0; then
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi


## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm -rf $SFTdir $outputv2_0 $outputv2_1 $outputv2_2 $outputv1_0 $outputv1_1 $outputv1_2 Fstats Fstats.log ${gridFile} ${outputv2_6}
fi
