#!/bin/bash

testDIR="./mfd_TEST"
  
oldcode=./makefakedata_test
newcodeDEFAULT=./lalapps_Makefakedata
compCode=./compareSFTs

# test if LAL_DATA_PATH has been set ... needed to locate ephemeris-files
if [ -z "$LAL_DATA_PATH" ]; then
    echo
    echo "Need environment-variable LAL_DATA_PATH to point to your ephemeris-directory (e.g. /usr/local/share/lal)"
    if [ -n "$LAL_PREFIX" ]; then
	echo "You have LAL_PREFIX set, I suggest setting 'LAL_DATA_PATH=\$LAL_PREFIX/share/lal'"
    fi
    echo
    exit 1
fi

if [ -z "$1" ]; then
    newcode=${newcodeDEFAULT}
else
    newcode="$1"
    echo "'$newcode'"
fi

#prepare test subdirectory
if [ ! -d "$testDIR" ]; then
    mkdir $testDIR
else
## cleanup: remove previous output-SFTs
    rm $testDIR/SFTtest_v2.00* || true
    rm $testDIR/SFTtest_v4.00* || true
fi

if [ ! -x "$oldcode" ]; then
    echo "Cannot run reference-code $oldcode!"
    exit -1;
elif [ ! -x "$newcode" ]; then
    echo "Cannot run test-code $newcode!"
    exit -1;
fi

# input parameters
## FIXED
ephemdir=${LAL_PREFIX}/share/lal
Tsft=1800
nTsft=20
timestamps=./testT8_1800
refTime=731210229
fmin=300.0
Band=10.0

## VARY
IFO=LLO
aPlus=1.5
aCross=0.7
psi=0.5
phi0=0.9
f0=300.2
alpha=1.7
delta=0.9
noiseDir="../"
noiseSFTs="$noiseDir/SFT.0000[0-9]"

dataTMP=In.data-test
oldCL="-i $dataTMP  -I $IFO -E $ephemdir -S $refTime -n ${testDIR}/SFTtest_v2" ## -D $noiseDir"
newCL="--Tsft=$Tsft --fmin=$fmin --Band=$Band --aPlus=$aPlus --aCross=$aCross --psi=$psi --phi0=$phi0 --f0=$f0 --latitude=$delta  --longitude=$alpha --detector=$IFO --outSFTbname=${testDIR}/SFTtest_v4 --timestampsFile=$timestamps --refTime=$refTime" ## -D$noiseSFTs -v1"

## produce In.data file for makefakedata_v2
echo "$Tsft	%Tsft_in_sec
$nTsft	%nTsft
$fmin   %first_SFT_frequency_in_Hz	
$Band	%SFT_freq_band_in_Hz
0.0	%sigma_(std_of_noise.When=0_only_signal_present)
$aPlus	%Aplus
$aCross	%Across
$psi	%psi
$phi0	%phi0
$f0	%f0
$delta	%latitude_in_radians
$alpha	%longitude_in_radians
0	%max_spin-down_param_order
$timestamps 	%name_of_time-stamps_file
" > In.data-test

echo "1) Testing isolated pulsar-signal without noise"
echo
echo "Running 'reference-code':"
echo "$oldcode $oldCL"

time $oldcode $oldCL

# remove temporary In.data - file
#rm $dataTMP
echo
echo "Running new Makefakedata code:"
echo "$newcode $newCL"
time $newcode $newCL

echo
echo "comparison of resulting SFTs:"

cmdline="$compCode -v -1 '${testDIR}/SFTtest_v2*' -2 '${testDIR}/SFTtest_v4*'"
echo ${cmdline}
eval ${cmdline}

