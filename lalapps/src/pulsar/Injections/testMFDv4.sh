#!/bin/bash

## allow 'make test' to work from builddir != srcdir
if [ -n "${srcdir}" ]; then
    builddir="./";
    msftdir="${builddir}../MakeSFTs/"
else
    srcdir=.
fi

mfdCODE="${builddir}lalapps_Makefakedata_v4"
cmpCODE="${builddir}lalapps_compareSFTs"
extractCODE="${msftdir}lalapps_ConvertToSFTv2"

testDIR1="./mfdv4_TEST1"
testDIR2="./mfdv4_TEST2"
testDIR3="./mfdv4_TEST3"


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

#prepare test subdirectory
if [ ! -d "$testDIR1" ]; then
    mkdir $testDIR1
else
## cleanup: remove previous output-SFTs
    rm -f $testDIR1/* || true
fi
if [ ! -d "$testDIR2" ]; then
    mkdir $testDIR2
else
## cleanup: remove previous output-SFTs
    rm -f $testDIR2/* || true
fi
if [ ! -d "$testDIR3" ]; then
    mkdir $testDIR3
else
## cleanup: remove previous output-SFTs
    rm -f $testDIR3/* || true
fi



tol="1e-4";	## tolerance on relative difference between SFTs in comparison
# input parameters
## FIXED
Tsft=1800
nTsft=20
timestamps="$srcdir/testT8_1800"
refTime=701210229

## excercise non-integer cycle gaps in heterodyned timeseries
fmin=299.1001
Band=9.9998
fmax=$(echo $fmin $Band | awk '{printf "%.7g", $1 + $2}');
fUpper=311

## VARY
IFO=H1
aPlus=1.5
aCross=0.7
psi=0.5
phi0=0.9
f0=300.2
alpha=1.7
delta=0.9

f1dot="-1.e-9"
f2dot="1e-14"


echo "------------------------------------------------------------"
echo " SIGNAL-ONLY - compare heterodyned SFTs with 'exact' ones"
echo "------------------------------------------------------------"

echo
echo "mfd_v4: producing SFTs via heterodyned timeseries (generationMode=0 [ALL_AT_ONCE] )..."
echo

mfdCL="--Tsft=$Tsft --fmin=$fmin --Band=$Band --aPlus=$aPlus --aCross=$aCross --psi=$psi --phi0=$phi0 --f0=$f0 --longitude=$alpha --latitude=$delta --detector=$IFO --timestampsFile=$timestamps --refTime=$refTime --f1dot=$f1dot --f2dot=$f2dot --generationMode=0 --outSFTbname=$testDIR1"
cmdline="$mfdCODE $mfdCL";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdCODE' ..."
    exit 1
fi


echo
echo "mfd_v4: producing SFTs via heterodyned timeseries (generationMode=1 [PER_SFT] )..."
echo

mfdCL="--Tsft=$Tsft --fmin=$fmin --Band=$Band --aPlus=$aPlus --aCross=$aCross --psi=$psi --phi0=$phi0 --f0=$f0 --longitude=$alpha --latitude=$delta --detector=$IFO --timestampsFile=$timestamps --refTime=$refTime --f1dot=$f1dot --f2dot=$f2dot --generationMode=1 --outSFTbname=$testDIR2"
cmdline="$mfdCODE $mfdCL";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdCODE' ..."
    exit 1
fi


echo
echo "mfd_v4: producing SFTs via 'exact' timeseries (non-heterodyned)..."
echo

mfdCL="--Tsft=$Tsft --fmin=0 --Band=$fUpper --aPlus=$aPlus --aCross=$aCross --psi=$psi --phi0=$phi0 --f0=$f0 --longitude=$alpha --latitude=$delta --detector=$IFO --timestampsFile=$timestamps --refTime=$refTime --f1dot=$f1dot --f2dot=$f2dot --generationMode=0 --outSFTbname=$testDIR3"
cmdline="$mfdCODE $mfdCL";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdCODE' ..."
    exit 1
fi


echo "... and extracting the relevant frequency band ..."
echo

## extract relevant frequency-band
extractCL="--inputSFTs='$testDIR3/*.sft' --fmin=$fmin --fmax=$fmax --outputDir=$testDIR3 --descriptionMisc=Band"
cmdline="$extractCODE $extractCL"
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$extractCODE' ..."
    exit 1
fi





echo
echo "comparison of resulting SFTs:"

cmdline="$cmpCODE -e $tol -1 '${testDIR1}/*.sft' -2 '${testDIR3}/*_Band*'"
echo ${cmdline}
if ! eval $cmdline; then
    echo "OUCH... SFTs differ by more than $tol. Something might be wrong..."
    exit 2
else
    echo "OK."
fi


echo
cmdline="$cmpCODE -e $tol -1 '${testDIR2}/*.sft' -2 '${testDIR3}/*_Band*'"
echo ${cmdline}
if ! eval $cmdline; then
    echo "OUCH... SFTs differ by more than $tol. Something might be wrong..."
    exit 2
else
    echo "OK."
fi


## clean up files [allow turning off via 'NOCLEANUP' environment variable
if [ -z "$NOCLEANUP" ]; then
    rm -rf $testDIR1 $testDIR2 $testDIR3
fi


