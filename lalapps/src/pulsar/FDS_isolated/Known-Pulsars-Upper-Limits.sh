#!/bin/bash

export LC_ALL=C

#------------------------------
#Some Constants
pi=3.14159265358979
pihalf=1.57079632679490
minuspihalf=-1.57079632679490
twopi=6.28318530717959

# Some codes which are needed for this run
randCode=makerandparam
mfd_code="lalapps_Makefakedata_v4"
cfsv2_code="lalapps_ComputeFStatistic_v2"

SFTdir="./Known-Pulsars-Upper-Limits_sfts"

maxiter=100;
echo "maxiter=$maxiter"
## ---------- param-ranges for MC
refTime=582370572;
fmin=411.01
Band=0.1
mfd_FreqBand=0.12;

## NOTE!!: use "0", NOT "0.0" !!!
f1min=-4.24E-16
f1max=-4.20E-16

## ---------- some derived input params
f0min=$(echo $fmin | awk '{printf "%.5f", $1 + 0.001}');
f0max=$(echo $fmin $Band | awk '{printf "%.5f", $1 + $2 - 0.001}');

#Tsft=1800;
#startTime=711595933
#refTime=711595933     ##$startTime
#duration=144000	      ## 27.7777 hours

IFO=H1

## --------------------------------------------------
## little helper function which converts a number from [0,1) to a given range [a,b)
## $1 input number from [0,1)
## $2 $3 output-range
## the result is stored in the global variable convToRangeOut
convToRangeOut="unset"
convToRange()
{
    if [ -z "$1" -o -z "$2" -o -z "$3" -o -n "$4" ]; then
	echo "Error convToRange() called with wrong number of arguments"
	exit 1
    fi

    local numin=$1
    local mini=$2
    local maxi=$3

    local numout=$(echo $numin $mini $maxi | awk '{printf "%g",($3-$2)*$1+$2}');

    ##'return' result in global variable
    convToRangeOut=$numout

} ## convToRange()

    h0=1.9e-25;

    Alpha=0.13289;
    Delta=0.08484;

iteration=1
#echo "psi          phi0        f0        f1dot        2F(CFS)" > resH1rand.txt
echo " Alpha         Delta         aPlus         aCross         cosiota         h0        psi         phi0        f0       f1dot" >> pulsar-411.06Hz.txt

while [ $iteration -le $maxiter ]; do
#echo
#echo "--------------------------------------"
echo "Iter = $iteration"
#echo "--------------------------------------"
## ----------------------------------------------------------------------
## generate random-values for the variable params:
## delta, alpha, IFO, aPlus, aCross, psi, phi0, f0

## generate deterministic random-seed:
    randseed=$((10 * $iteration + $iteration))
    randvals=`./$randCode -s $randseed -n 9`

     #---------- psi
    randval=$(echo $randvals | awk '{print $6}');
    convToRange $randval 0 $pi
    psi=$convToRangeOut
    #---------- phi0
    randval=$(echo $randvals | awk '{print $7}');
    convToRange $randval 0 $twopi
    phi0=$convToRangeOut
    #---------- f0
    randval=$(echo $randvals | awk '{print $8}');
    convToRange $randval $f0min $f0max
    f0=$convToRangeOut
    #---------- f1dot
    randval=$(echo $randvals | awk '{print $9}');
    convToRange $randval $f1min $f1max
    f1dot=$convToRangeOut
    #---------- f2dot
    randval=$(echo $randvals | awk '{print $10}');
    convToRange $randval -1 1
    cosiota=$convToRangeOut

    aPlus=$(echo $h0 $cosiota | awk '{printf "%1.5e", 0.5*$1*(1 + $2 * $2)}');
    aCross=$(echo $h0 $cosiota | awk '{printf "%1.5e", $1*$2}');

#--------------------------------------------------
## test starts here
##--------------------------------------------------

#echo
#echo "----------------------------------------------------------------------"
#echo " STEP 1: Generate Fake Signal"
#echo "----------------------------------------------------------------------"
#echo
#echo
if [ ! -d "$SFTdir" ]; then
    mkdir $SFTdir;
else
    rm -f $SFTdir/H-1_H1_1800SFT_mfdv4-81*;
    rm -f $SFTdir/H-1_H1_1800SFT_mfdv4-82*;
    rm -f $SFTdir/H-1_H1_1800SFT_mfdv4-83*;

fi

mfd_fmin=$(echo $f0 $mfd_FreqBand | awk '{printf "%g", $1 - $2/2.0}');

#echo "----------------------------------------------------------------------"
# concatenate this with the mfd-specific switches:
mfd_CL="--detector=$IFO --latitude=$Delta  --longitude=$Alpha --aPlus=$aPlus --aCross=$aCross --psi=$psi --phi0=$phi0 --fmin=$mfd_fmin --Band=$mfd_FreqBand --f0=$f0 --outSFTbname=$SFTdir/ --f1dot=$f1dot --timestampsFile='./TimeStamps.txt' --generationMode=1 --refTime=$refTime --ephemYear=05-09 --noiseSFTs='/local_data/gholami/LAL-sources/lalapps-CVS/src/pulsar/FDS_isolated/workdir/Pulsars/Known-Pulsars/S5-SFTs-410.9-411.2Hz/H-1_H1*.sft'"
############################################################
cmdline="$mfd_code $mfd_CL";
#echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

#echo "----------------------------------------------------------------------"

#FreqBand=$(echo $mfd_FreqBand | awk '{printf "%g", ($1-0.05*$1) / 2.0}');

CFSv2_CL="--Freq=$f0 --f1dot=$f1dot --Alpha=$Alpha --Delta=$Delta --ephemYear=05-09 --minStartTime=818881041 --maxEndTime=835359537 --TwoFthreshold=0 --refTime=$refTime --DataFiles='$SFTdir/*.sft' --outputFstat='Known-Pulsars-Upper-Limit-411.06Hz-h0=$h0-$maxiter Trials.txt' -H -A"
############################################################
cmdline="$cfsv2_code $CFSv2_CL";
#echo $cmdline;

if ! eval $cmdline; then
    echo "Error.. something failed when running '$cfsv2_code' ..."
    exit 1
fi

echo " $Alpha    $Delta    $aPlus    $aCross    $cosiota     $h0    $psi     $phi0    $f0   $f1dot" >> pulsar-411.06Hz.txt

    iteration=$(($iteration + 1));
done
## end of while-loop
