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
cfs_code_orig="${builddir}lalapps_ComputeFStatistic_v2_orig"
cfs_code_test="${builddir}lalapps_ComputeFStatistic_v2"
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

Ftolerance=0.01
# ---------- fixed parameter of our test-signal
Tsft=1800;
startTime=711595934
refTime=701595833  ## $startTime
duration=144000		## 40 hours

mfd_FreqBand=5.0;

Alpha=2.0
Delta=-0.5

h0=0.1
cosi=-0.3

psi=0.6
phi0=1.5

Freq=1000.12345
mfd_fmin=$(echo $Freq $mfd_FreqBand | awk '{printf "%g", $1 - $2 / 2.0}');

f1dot=-1e-7;
f2dot=1e-14

noiseSqrtSh=1

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
mfd_CL="${saf_CL} --fmin=$mfd_fmin --Band=$mfd_FreqBand --Freq=$Freq --outSFTbname=$SFTdir --f1dot=$f1dot --f2dot=$f2dot --refTime=$refTime"
if [ "$haveNoise" = true ]; then
    mfd_CL="$mfd_CL --noiseSqrtSh=$sqrtSh"; ## --randSeed=1";
fi

cmdline="$mfd_code $mfd_CL";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo "STEP 2: run original CFS_v2 with perfect match"
echo "----------------------------------------------------------------------"
echo
outfile_orig="Fstat_v2_orig.dat";
## common cmdline-options for v1 and v2
cfs_CL="--IFO=$IFO --Freq=$Freq --Alpha=$Alpha --Delta=$Delta --f1dot=$f1dot --f2dot=$f2dot --DataFiles='$SFTdir/*.sft' --refTime=$refTime --TwoFthreshold=0"
if [ "$haveNoise" = false ]; then
    cfs_CL="$cfs_CL --SignalOnly"
fi

cmdline="$cfs_code_orig $cfs_CL  --outputFstat=$outfile_orig";
echo $cmdline;

if ! eval $cmdline; then
    echo "Error.. something failed when running '$cfs_code_orig' ..."
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo " STEP 3: run REAL4 CFS_v2 with perfect match"
echo "----------------------------------------------------------------------"
echo
outfile_test="Fstat_v2_test.dat";
cmdline="$cfs_code_test $cfs_CL --outputFstat=$outfile_test";
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cfs_code_test' ..."
    exit 1;
fi

echo
echo "----------------------------------------"
echo " STEP 4: Comparing results: "
echo "----------------------------------------"
echo
echo "----- CFS_v2_orig: "
cat $outfile_orig | sed -e"/^%.*/{d}"
echo "----- CFS_v2_test: "
cat $outfile_test | sed -e"/^%.*/{d}"

echo
cmdline="$cmp_code -1 ./$outfile_orig -2 ./$outfile_test --clusterFiles=0 --Ftolerance=$Ftolerance"
echo $cmdline
if ! eval $cmdline; then
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi

## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm -rf $SFTdir $outfile_orig $outfile_test Fstats Fstats.log
fi
