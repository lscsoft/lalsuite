#!/bin/bash

#prepare test subdirectory
if [ ! -d testMakefakedata ]; then
    mkdir testMakefakedata && cd testMakefakedata;
else
## cleanup: remove previous output-SFTs
    cd testMakefakedata
    rm SFTtest_v2.00* || 1
    rm SFTtest_v4.00* || 1
fi
  
oldcode=../makefakedata_test
newcode=../lalapps_Makefakedata

if [ ! -x $oldcode ]; then
    echo "Cannot run reference-code $oldcode!"
    exit -1;
elif [ ! -x $newcode ]; then
    echo "Cannot run test-code $newcode!"
    exit -1;
fi

# signal parameters
IFO=LLO
ephemdir=../ephems
startTime=731210229
refTime=$startTime
Tsft=1800
nTsft=10
fmin=1008.0
Band=3.0
aPlus=1.5
aCross=0.7
psi=0.5
phi0=0.9
f0=1009.0
alpha=1.7
delta=0.9

dataTMP=In.data-test
oldCL="-i $dataTMP  -I $IFO -E $ephemdir -G $startTime -S $refTime -n SFTtest_v2"
newCL="--Tsft=$Tsft --nTsft=$nTsft --fmin=$fmin --Band=$Band --aPlus=$aPlus --aCross=$aCross --psi=$psi --phi0=$phi0 --f0=$f0 --latitude=$delta  --longitude=$alpha --detector=$IFO --ephemDir=$ephemdir --outSFTbname=SFTtest_v4 --startTime=$startTime --refTime=$refTime"

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
T8	%name_of_time-stamps_file
" > In.data-test

echo "1) Testing isolated pulsar-signal without noise"
echo "Running 'reference-code' $oldcode":
$oldcode $oldCL &> /dev/null

# remove temporary In.data - file
#rm $dataTMP

echo "Running new code $newcode:"
$newcode $newCL
echo "comparison:"
../compareSFTs -1 "./SFTtest_v2*" -2 "./SFTtest_v4*"


cd ..
