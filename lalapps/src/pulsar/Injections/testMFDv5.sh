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
fmin=299.0
Band=10
fmax=$(echo $fmin $Band | LC_ALL=C awk '{printf "%.7g", $1 + $2}');

IFO1=H1
IFO2=L1
sqrtSn1=0.5;
sqrtSn2=1.2;	## for comparison with 2 calls to mfdv4 and fixed see, the 2nd IFO noise must be 0

timestamps1="${testDIR}/H1-timestamps.dat"
timestamps2="${testDIR}/L1-timestamps.dat"

cp ${timestamps} ${timestamps1}
echo "701210229 0" > ${timestamps2}
echo "701230229 0" >> ${timestamps2}
echo "701240229 0" >> ${timestamps2}

## ---------- signal parameters ----------
## ----- signal 1
s1_h0=0.73
s1_cosi=0.1
s1_psi=0.5
s1_phi0=0.9
s1_refTime=701210229
s1_Alpha=1.7
s1_Delta=0.9
s1_Freq=299.12
s1_f1dot=0
s1_f2dot=0
## ----- signal 2
s2_h0=2.5
s2_cosi=-0.5
s2_psi=1.2
s2_phi0=1.5
s2_refTime=711210229
s2_Alpha=3.7
s2_Delta=-0.5
s2_Freq=300.12
s2_f1dot=-1e-10
s2_f2dot=0
## ----- signal 3
s3_h0=3.1
s3_cosi=0.5
s3_psi=-1.2
s3_phi0=2.5
s3_refTime=721210229
s3_Alpha=0.5
s3_Delta=1.2
s3_Freq=300.00
s3_f1dot=-1e-9
s3_f2dot=-2e-19
# --------------------

injString="Alpha=${s1_Alpha};Delta=${s1_Delta};refTime=${s1_refTime};Freq=${s1_Freq};f1dot=${s1_f1dot};f2dot=${s1_f2dot};h0=${s1_h0};cosi=${s1_cosi};psi=${s1_psi};phi0=${s1_phi0};"

## ---------- signal file 1 ----------
injFile1=${testDIR}/injectionS1.dat
echo "[Pulsar 1]" >> ${injFile1}
echo "Alpha = ${s1_Alpha}" >> ${injFile1}
echo "Delta = ${s1_Delta}" >> ${injFile1}
echo "refTime = ${s1_refTime}" >> ${injFile1}
echo "Freq = ${s1_Freq}" >> ${injFile1}
echo "f1dot = ${s1_f1dot}" >> ${injFile1}
echo "f2dot = ${s1_f2dot}" >> ${injFile1}
echo "h0 = ${s1_h0}" >> ${injFile1}
echo "cosi = ${s1_cosi}" >> ${injFile1}
echo "psi = ${s1_psi}" >> ${injFile1}
echo "phi0 = ${s1_phi0}" >> ${injFile1}
echo >> ${injFile1}
## ---------- signal file 2 ----------
injFile2=${testDIR}/injectionS2.dat
echo "Alpha = ${s2_Alpha}" >> ${injFile2}
echo "Delta = ${s2_Delta}" >> ${injFile2}
echo "refTime = ${s2_refTime}" >> ${injFile2}
echo "Freq = ${s2_Freq}" >> ${injFile2}
echo "f1dot = ${s2_f1dot}" >> ${injFile2}
echo "f2dot = ${s2_f2dot}" >> ${injFile2}
echo "h0 = ${s2_h0}" >> ${injFile2}
echo "cosi = ${s2_cosi}" >> ${injFile2}
echo "psi = ${s2_psi}" >> ${injFile2}
echo "phi0 = ${s2_phi0}" >> ${injFile2}
echo >> ${injFile2}
## ---------- add section for Pulsar 3 into signal-file 2 ----------
echo "[Pulsar 3]" >> ${injFile2}
echo "Alpha = ${s3_Alpha}" >> ${injFile2}
echo "Delta = ${s3_Delta}" >> ${injFile2}
echo "refTime = ${s3_refTime}" >> ${injFile2}
echo "Freq = ${s3_Freq}" >> ${injFile2}
echo "f1dot = ${s3_f1dot}" >> ${injFile2}
echo "f2dot = ${s3_f2dot}" >> ${injFile2}
echo "h0 = ${s3_h0}" >> ${injFile2}
echo "cosi = ${s3_cosi}" >> ${injFile2}
echo "psi = ${s3_psi}" >> ${injFile2}
echo "phi0 = ${s3_phi0}" >> ${injFile2}
echo >> ${injFile2}


## ---------- output parameters ----------
sftsv4_1=${testDIR}/${IFO1}-sftsv4.sft
sftsv4_2=${testDIR}/${IFO2}-sftsv4.sft
sftsv5_1=${testDIR}/H-*_mfdv5*.sft
sftsv5_2=${testDIR}/L-*_mfdv5*.sft

## ----------
## produce SFTs for 2 detectors, containing Gaussian noise + N signals, compare between mfdv4 and mfdv5
## ----------

echo
echo "========== MFDv4 =========="
echo
mfdv4_CL="$mfdv4_CODE ${mfdv4_extra} --fmin=$fmin --Band=$Band --generationMode=0  -v${debug} --outSingleSFT"
echo "---------- mfdv4: inject first signal ----------"
sig1="--refTime=${s1_refTime} --h0=${s1_h0} --cosi=${s1_cosi} --psi=${s1_psi} --phi0=${s1_phi0} --Freq=${s1_Freq} --Alpha=${s1_Alpha} --Delta=${s1_Delta} --f1dot=${s1_f1dot} --f2dot=${s1_f2dot}"
##----- first IFO
out_IFO1="--IFO=${IFO1} --timestampsFile=${timestamps1}  --outSFTbname=${sftsv4_1} --noiseSqrtSh=${sqrtSn1} --randSeed=1"
cmdline="$mfdv4_CL ${sig1} ${out_IFO1}"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdv4_CODE' ..."
    exit 1
fi
##----- second IFO
out_IFO2="--IFO=${IFO2} --timestampsFile=${timestamps2}  --outSFTbname=${sftsv4_2} --noiseSqrtSh=${sqrtSn2} --randSeed=2"
cmdline="$mfdv4_CL  ${sig1} ${out_IFO2}"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdv4_CODE' ..."
    exit 1
fi
echo "---------- mfdv4: inject second signal on top ----------"
sig2="--refTime=${s2_refTime} --h0=${s2_h0} --cosi=${s2_cosi} --psi=${s2_psi} --phi0=${s2_phi0} --Freq=${s2_Freq} --Alpha=${s2_Alpha} --Delta=${s2_Delta} --f1dot=${s2_f1dot} --f2dot=${s2_f2dot}"
##----- first IFO
out_IFO1="--IFO=${IFO1} --noiseSFTs=${sftsv4_1} --window=None --outSFTbname=${sftsv4_1}"
cmdline="$mfdv4_CL ${sig2} ${out_IFO1}"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdv4_CODE' ..."
    exit 1
fi
##----- second IFO
out_IFO2="--IFO=${IFO2} --noiseSFTs=${sftsv4_2} --window=None --outSFTbname=${sftsv4_2}"
cmdline="$mfdv4_CL ${sig2} ${out_IFO2}"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdv4_CODE' ..."
    exit 1
fi
echo "---------- mfdv4: inject third signal on top ----------"
sig2="--refTime=${s3_refTime} --h0=${s3_h0} --cosi=${s3_cosi} --psi=${s3_psi} --phi0=${s3_phi0} --Freq=${s3_Freq} --Alpha=${s3_Alpha} --Delta=${s3_Delta} --f1dot=${s3_f1dot} --f2dot=${s3_f2dot}"
##----- first IFO
out_IFO1="--IFO=${IFO1} --noiseSFTs=${sftsv4_1} --window=None --outSFTbname=${sftsv4_1}"
cmdline="$mfdv4_CL ${sig2} ${out_IFO1}"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdv4_CODE' ..."
    exit 1
fi
##----- second IFO
out_IFO2="--IFO=${IFO2} --noiseSFTs=${sftsv4_2} --window=None --outSFTbname=${sftsv4_2}"
cmdline="$mfdv4_CL ${sig2} ${out_IFO2}"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdv4_CODE' ..."
    exit 1
fi


echo
echo "========== MFDv5 =========="
echo
mfdv5_CL="$mfdv5_CODE ${mfdv5_extra} --outSingleSFT --outSFTdir=${testDIR} --fmin=$fmin --Band=$Band -v${debug}"

echo "----- single multi-IFO, multi-signal call"
outIFOs="--IFOs=${IFO1},${IFO2} --timestampsFiles=${timestamps1},${timestamps2} --sqrtSX=${sqrtSn1},${sqrtSn2} --randSeed=1"
sig13="--injectionSources='@${injFile1};${injFile2}'"
cmdline="$mfdv5_CL ${outIFOs} ${sig13}"
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdv5_CODE' ..."
    exit 1
fi

echo "----- and again the same, using different input methods"
sig1="--injectionSources='${injString}'"
cmdline="$mfdv5_CL ${outIFOs} ${sig1} --outMiscField='mfdv5_try2'"
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdv5_CODE' ..."
    exit 1
fi


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
