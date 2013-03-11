#!/bin/bash

## allow 'make test' to work from builddir != srcdir
if [ -z "${srcdir}" ]; then
    srcdir=`dirname $0`
fi

builddir="./";
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
    mfdv4_extra="-E ${LALPULSAR_DATADIR}"
    mfdv5_extra="-E ${LALPULSAR_DATADIR}"
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


tol=1e-3;	## tolerance on relative difference between SFTs in comparison
# input parameters
## ---------- data parameters ----------
Tsft=1800
nTsft=20
timestamps=${srcdir}/testT8_1800

## excercise non-integer cycle gaps in heterodyned timeseries
fmin=299.1001
Band=0.12345
fmax=$(echo $fmin $Band | LC_ALL=C awk '{printf "%.7g", $1 + $2}');

IFO1=H1
IFO2=L1
IFOs=${IFO1},${IFO2}

sqrtSn1=0;
sqrtSn2=0;	## for comparison with 2 calls to mfdv4 and fixed see, the 2nd IFO noise must be 0
sqrtSnX=${sqrtSn1},${sqrtSn2}

timestamps1="${testDIR}/H1-timestamps.dat"
timestamps2="${testDIR}/L1-timestamps.dat"
timestampsFiles=${timestamps1},${timestamps2}

cp ${timestamps} ${timestamps1}
echo "701210229 0" > ${timestamps2}
echo "701230229 0" >> ${timestamps2}
echo "701240229 0" >> ${timestamps2}

## ---------- signal parameters ----------
h0=0.73
cosi=0.1
psi=0.5
phi0=0.9

refTime=701210229
Alpha=1.7
Delta=0.9
Freq=299.12
f1dot=0
f2dot=0

injFile=${testDIR}/injectionSources.dat
echo "Alpha = ${Alpha}" > ${injFile}
echo "Delta = ${Delta}" >> ${injFile}
echo "refTime = ${refTime}" >> ${injFile}
echo "Freq = ${Freq}" >> ${injFile}
echo "f1dot = ${f1dot}" >> ${injFile}
echo "f2dot = ${f2dot}" >> ${injFile}

echo "h0 = ${h0}" >> ${injFile}
echo "cosi = ${cosi}" >> ${injFile}
echo "psi = ${psi}" >> ${injFile}
echo "phi0 = ${phi0}" >> ${injFile}

## ---------- output parameters ----------
sftsv4_1=${testDIR}/${IFO1}-sftsv4.sft
sftsv4_2=${testDIR}/${IFO2}-sftsv4.sft
sftsv5_1=${testDIR}/H-*_mfdv5*.sft
sftsv5_2=${testDIR}/L-*_mfdv5*.sft

echo "------------------------------------------------------------"
echo " SIGNAL-ONLY - compare SFTs between mfd_v4 and mfd_v5"
echo "------------------------------------------------------------"

echo
echo "----- mfd_v4: producing SFTs via (generationMode=1 [PER_SFT] ):"
echo
##----- first IFO
mfdv4_CL="$mfdv4_CODE ${mfdv4_extra} --Tsft=$Tsft --fmin=$fmin --Band=$Band --h0=$h0 --cosi=$cosi --psi=$psi --phi0=$phi0 --Freq=${Freq} --Alpha=$Alpha --Delta=$Delta --IFO=$IFO1 --timestampsFile=$timestamps1 --refTime=$refTime --f1dot=$f1dot --f2dot=$f2dot --generationMode=1 --noiseSqrtSh=${sqrtSn1} --randSeed=1 -v${debug} --outSingleSFT --outSFTbname=${sftsv4_1}"
echo $mfdv4_CL;
if ! eval $mfdv4_CL; then
    echo "Error.. something failed when running '$mfdv4_CODE' ..."
    exit 1
fi
echo "ok."

##----- second IFO
mfdv4_CL="$mfdv4_CODE ${mfdv4_extra} --Tsft=$Tsft --fmin=$fmin --Band=$Band --h0=$h0 --cosi=$cosi --psi=$psi --phi0=$phi0 --Freq=${Freq} --Alpha=$Alpha --Delta=$Delta --IFO=$IFO2 --timestampsFile=$timestamps2 --refTime=$refTime --f1dot=$f1dot --f2dot=$f2dot --generationMode=1 --noiseSqrtSh=${sqrtSn2} --randSeed=1 -v${debug} --outSingleSFT --outSFTbname=${sftsv4_2}"
echo $mfdv4_CL;
if ! eval $mfdv4_CL; then
    echo "Error.. something failed when running '$mfdv4_CODE' ..."
    exit 1
fi
echo "ok."


echo
echo "----- mfd_v5: producing SFTs:"
echo
## ----- multi-IFO call
mfdv5_CL="$mfdv5_CODE ${mfdv5_extra} --outSingleSFT --outSFTdir=${testDIR} --Tsft=$Tsft --fmin=$fmin --Band=$Band --IFOs=${IFOs} --timestampsFiles=${timestampsFiles} --sqrtSX=${sqrtSnX} --randSeed=1 --injectionSources=@${injFile} -v${debug}"
cmdline=" $mfdv5_CL "
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdv5_CODE' ..."
    exit 1
fi
echo "ok."

echo
echo "--------------------------------------------------"
echo "Comparison of resulting (concatenated) SFTs:"
echo "--------------------------------------------------"

cmdline="$cmp_CODE -e ${tol} -1 ${sftsv4_1} -2 '${sftsv5_1}' -d${debug}"
echo ${cmdline}
if ! eval $cmdline; then
    echo "Failed. SFTs produced by makefakedata_v4 and makefakedata_v5 differ by more than ${tol}!"
    exit 2
else
    echo "OK."
fi

cmdline="$cmp_CODE -e ${tol} -1 ${sftsv4_2} -2 '${sftsv5_2}' -d${debug}"
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
