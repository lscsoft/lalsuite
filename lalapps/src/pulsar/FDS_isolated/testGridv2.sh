#!/bin/sh

## set LAL debug level
echo "Setting LAL_DEBUG_LEVEL=${LAL_DEBUG_LEVEL:-msglvl1,memdbg}"
export LAL_DEBUG_LEVEL

## allow 'make test' to work from builddir != srcdir
if [ -z "${srcdir}" ]; then
    srcdir=`dirname $0`
fi

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
cfsv2_code="${builddir}lalapps_ComputeFStatistic_v2"
cmp_code="${builddir}lalapps_compareFstats"

SFTdir="./testGridv2_sfts"

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
echo " STEP 2: run CFS_v2 for various grid-types"
echo "----------------------------------------------------------------------"
echo
gridType0=" --gridType=0";	## flat grid
gridType1=" --gridType=1";	## isotropic
gridType2=" --gridType=2 --metricType=1 --metricMismatch=0.1";	## metric grid
output_0="./testGridv2_grid0.dat";
output_1="./testGridv2_grid1.dat";
output_2="./testGridv2_grid2.dat";
output_6="./testGridv2_grid6.dat";

sky_CL="--Alpha=$Alpha --AlphaBand=$AlphaBand --dAlpha=$dAlpha --Delta=$Delta --DeltaBand=$DeltaBand --dDelta=$dDelta"
spin_CL="--Freq=$Freq --FreqBand=$FreqBand --dFreq=$dFreq --f1dot=$f1dot --f1dotBand=$f1dotBand --df1dot=$df1dot"
cfs_CL="--IFO=$IFO --DataFiles='$SFTdir/testSFT*'"
if [ "$haveNoise" = false ]; then
    cfs_CL="$cfs_CL --SignalOnly"
fi

cmdline="$cfsv2_code $cfs_CL --TwoFthreshold=0 --Dterms=16 $extra_args";

## ----- grid=0
echo "CFSv2 using gridType=0:"
cmd0="$cmdline ${sky_CL} ${spin_CL} $gridType0 --outputFstat=$output_0";
echo $cmd0
if ! eval $cmd0; then
    echo "Error.. something failed when running '$cmd0' ..."
    exit 1
fi

## ----- grid=1
echo "CFSv2 using gridType=1:"
cmd0="$cmdline ${sky_CL} ${spin_CL} $gridType1 --outputFstat=$output_1";
echo $cmd0
if ! eval $cmd0; then
    echo "Error.. something failed when running '$cmd0' ..."
    exit 1
fi

## ----- grid=2
echo "CFSv2 using gridType=2:"
cmd0="$cmdline ${sky_CL} ${spin_CL} $gridType2 --outputFstat=$output_2";
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
grid_line=$(sed '/^%.*/d' ${output_2} | awk "$awk_extract6")


cmd0="$cmdline --gridType=6 --gridFile=./${gridFile} --outputFstat=$output_6";
echo $cmd0
if ! eval $cmd0; then
    echo "Error.. something failed when running '$cmd2' ..."
    exit 1
fi


echo
echo "----------------------------------------"
echo " STEP 3: Compare to reference results: "
echo "----------------------------------------"
echo

## ----- grid=0
echo "Comparing gridType=0:"
cmd0="$cmp_code -1 ${output_0} -2 ${srcdir}/${output_0}.ref.gz --clusterFiles=0 --Ftolerance=0.1";
echo $cmd0
if ! eval $cmd0; then
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi

## ----- grid=1
echo "Comparing gridType=1:"
cmd0="$cmp_code -1 ${output_1} -2 ${srcdir}/${output_1}.ref.gz --clusterFiles=0 --Ftolerance=0.1";
echo $cmd0
if ! eval $cmd0; then
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi


## ----- grid=2
echo "Comparing gridType=2:"
cmd0="$cmp_code -1 ${output_2} -2 ${srcdir}/${output_2}.ref.gz --clusterFiles=0 --Ftolerance=0.1";
echo $cmd0
if ! eval $cmd0; then
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi

## ----- grid=6
echo "Comparing gridType=6:"
cmd0="$cmp_code -1 ${output_6} -2 ${srcdir}/${output_6}.ref.gz --clusterFiles=0 --Ftolerance=0.01";
echo $cmd0
if ! eval $cmd0; then
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi


## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm -rf $SFTdir $output_0 $output_1 $output_2 $output_6 Fstats Fstats.log $gridFile
fi
