#!/bin/bash

## allow 'make test' to work from builddir != srcdir
if [ -z "${srcdir}" ]; then
    srcdir=`dirname $0`
fi

builddir="./";
msftdir="${builddir}../MakeSFTs/"

mfdv4_CODE="${builddir}lalapps_Makefakedata_v4"
mfdv5_CODE="${builddir}lalapps_Makefakedata_v5"
cmp_CODE="${builddir}lalapps_compareSFTs"

testDIR="./mfdv5_TEST"


## ----- user-controlled level of debug-output detail
if [ -n "$DEBUG" ]; then
    debug=${DEBUG}
else
    debug=0	## default=quiet
fi

if [ -n "${LALPULSAR_DATADIR}" ]; then
    mfdv4_CODE="${mfdv4_CODE} -E ${LALPULSAR_DATADIR}"
    mfdv5_CODE="${mfdv5_CODE} -E ${LALPULSAR_DATADIR}"
else
    echo
    echo "Need environment-variable LALPULSAR_DATADIR to be set to"
    echo "your ephemeris-directory (e.g. /usr/local/share/lalpulsar)"
    echo "This might indicate an incomplete LAL+LALPULSAR installation"
    echo
    exit 1
fi

## cleanup: remove any previous output-SFTs
rm -rf ${testDIR} || true
#prepare test subdirectory
mkdir -p $testDIR


tol=1e-10;	## tolerance on relative difference between SFTs in comparison
# input parameters
## FIXED
Tsft=1800
nTsft=20
timestamps=${srcdir}/testT8_1800
refTime=701210229

## excercise non-integer cycle gaps in heterodyned timeseries
fmin=299.1001
Band=9.9998
fmax=$(echo $fmin $Band | LC_ALL=C awk '{printf "%.7g", $1 + $2}');
fUpper=311

## VARY
IFO=H1
h0=0.73
sqrtSn=1;
cosi=0.1
psi=0.5
phi0=0.9
Freq=300.2
alpha=1.7
delta=0.9

f1dot=-1.e-9
f2dot=1e-14


sftsv4=${testDIR}/sftsv4.sft
sftsv5=${testDIR}/sftsv5.sft

echo "------------------------------------------------------------"
echo " SIGNAL-ONLY - compare SFTs between mfd_v4 and mfd_v5"
echo "------------------------------------------------------------"

echo
echo "----- mfd_v4: producing SFTs via (generationMode=0 [ALL_AT_ONCE] ):"
echo

mfdv4_CL="--Tsft=$Tsft --fmin=$fmin --Band=$Band --h0=$h0 --cosi=$cosi --psi=$psi --phi0=$phi0 --Freq=${Freq} --Alpha=$alpha --Delta=$delta --IFO=$IFO --timestampsFile=$timestamps --refTime=$refTime --f1dot=$f1dot --f2dot=$f2dot --generationMode=0 --noiseSqrtSh=${sqrtSn} --randSeed=1 -v${debug}"
cmdline="$mfdv4_CODE $mfdv4_CL --outSingleSFT --outSFTbname=${sftsv4}";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdCODE' ..."
    exit 1
fi
echo "ok."

echo
echo "----- mfd_v5: producing SFTs:"
echo

mfdv5_CL="--Tsft=$Tsft --fmin=$fmin --Band=$Band --h0=$h0 --cosi=$cosi --psi=$psi --phi0=$phi0 --Freq=${Freq} --Alpha=$alpha --Delta=$delta --IFO=$IFO --timestampsFile=$timestamps --refTime=$refTime --f1dot=$f1dot --f2dot=$f2dot --noiseSqrtSh=${sqrtSn} --randSeed=1 -v${debug}"
cmdline="$mfdv5_CODE $mfdv5_CL --outSingleSFT --outSFTbname=${sftsv5}";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdCODE' ..."
    exit 1
fi
echo "ok."

echo
echo "--------------------------------------------------"
echo "Comparison of resulting (concatenated) SFTs:"
echo "--------------------------------------------------"

cmdline="$cmp_CODE -e ${tol} -1 ${sftsv4} -2 ${sftsv5} -d${debug}"
echo ${cmdline}
if ! eval $cmdline; then
    echo "Failed. SFTs produced by makefakedata_v4 and makefakedata_v5 differ by more than ${tol}!"
    exit 2
else
    echo "OK."
fi

## clean up files [allow turning off via 'NOCLEANUP' environment variable
if [ -z "$NOCLEANUP" ]; then
    rm -rf ${testDIR}
fi
