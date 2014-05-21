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
saf_code="lalapps_SemiAnalyticF"
mfd_code="lalapps_Makefakedata_v4"
cfs_code="lalapps_ComputeFStatistic"
cfsv2_code="lalapps_ComputeFStatistic_v2"
cmp_code="lalapps_compareFstats"
pfs_code="lalapps_PredictFStat"

SFTdir="./testRandom_sfts"

maxiter=1;
echo "maxiter=$maxiter"
## ---------- param-ranges for MC
fmin=50.0
Band=950.0
## NOTE!!: use "0", NOT "0.0" !!!
f1min=-1e-7
f1max=1e-7

## ---------- some derived input params
f0min=$(echo $fmin | awk '{printf "%.8f", $1 + 1}');
f0max=$(echo $fmin $Band | awk '{printf "%.8f", $1 + $2 - 1}');

Tsft=1800;
startTime=711595933
refTime=711595933     ##$startTime
duration=144000	      ## 27.7777 hours

mfd_FreqBand=1.0;

aplusmin=1e-21
aplusmax=1e-20
acrossmin=0.5e-21
acrossmax=0.5e-20

noiseSqrtSh=3e-20

if [ "$noiseSqrtSh" != 0 ]; then
    sqrtSh=$noiseSqrtSh
    whatNoise=
else
    sqrtSh=1;
    whatNoise="--SignalOnly";
fi

IFO1=H1
IFO2=H2
IFO3=L1
IFO4=G1

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


iteration=1
echo " Alpha        Delta        aPlus        aCross         psi          phi0        f0        f1dot        2F(PFS)     2F(CFS)" > resH1rand.txt
echo " Alpha        Delta        aPlus        aCross         psi          phi0        f0        f1dot        2F(PFS)     2F(CFS)" > resH2rand.txt
echo " Alpha        Delta        aPlus        aCross         psi          phi0        f0        f1dot        2F(PFS)     2F(CFS)" > resL1rand.txt
echo " Alpha        Delta        aPlus        aCross         psi          phi0        f0        f1dot        2F(PFS)     2F(CFS)" > resG1rand.txt
echo " Alpha        Delta        aPlus        aCross         psi          phi0        f0        f1dot        2F(PFS)     2F(CFS)" > resH1H2rand.txt
echo " Alpha        Delta        aPlus        aCross         psi          phi0        f0        f1dot        2F(PFS)     2F(CFS)" > resH1H2L1G1rand.txt


while [ $iteration -le $maxiter ]; do
#echo
#echo "--------------------------------------"
echo "Iter = $iteration"
#echo "--------------------------------------"
## ----------------------------------------------------------------------
## generate random-values for the variable params:
## delta, alpha, IFO, aPlus, aCross, psi, phi0, f0

## generate deterministic random-seed:
    randseed=$(( 1000 * $iteration + $iteration ))
    randvals=`$randCode -s $randseed -n 9`

    #---------- delta
    randval=$(echo $randvals | awk '{print $1}');
    Delta=$randval
    #---------- alpha
    randval=$(echo $randvals | awk '{print $2}');
    convToRange $randval 0 $twopi
    Alpha=$convToRangeOut
    #---------- aPlus
    randval=$(echo $randvals | awk '{print $4}');
    convToRange $randval $aplusmin $aplusmax
    aPlus=$convToRangeOut
    #---------- aCross
    randval=$(echo $randvals | awk '{print $5}');
    convToRange $randval $acrossmin $acrossmax
    aCross=$convToRangeOut
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

##--------------------------------------------------
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
    rm -f $SFTdir/*;
fi

mfd_fmin=$(echo $f0 $mfd_FreqBand | awk '{printf "%g", $1 - $2 / 2.0}');

# this part of the command-line is compatible with SemiAnalyticF:
saf_CL="--latitude=$Delta  --longitude=$Alpha --Tsft=$Tsft --startTime=$startTime --duration=$duration --aPlus=$aPlus --aCross=$aCross --psi=$psi --phi0=$phi0"
# concatenate this with the mfd-specific switches:
mfd_CL="${saf_CL} --fmin=$mfd_fmin --Band=$mfd_FreqBand --f0=$f0 --outSFTbname=$SFTdir/ --f1dot=$f1dot  --refTime=$refTime --noiseSqrtSh=$sqrtSh"
############################################################
cmdline="$mfd_code $mfd_CL --detector=$IFO1";
#echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi
############################################################
cmdline="$mfd_code $mfd_CL --detector=$IFO2";
#echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi
############################################################
cmdline="$mfd_code $mfd_CL --detector=$IFO3";
#echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi
############################################################
cmdline="$mfd_code $mfd_CL --detector=$IFO4";
#echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

#echo
#echo "----------------------------------------------------------------------"
#echo " STEP 2: Calculate F-Statistic Semi-Analytically"
#echo "----------------------------------------------------------------------"
# this part of the command-line is compatible with PredictFStat:
FreqBand=$(echo $mfd_FreqBand | awk '{printf "%g", ($1-0.05*$1) / 2.0}');
pfs_CL="--aPlus=$aPlus --aCross=$aCross --psi=$psi --Freq=$f0 --Delta=$Delta --Alpha=$Alpha"

############################################################
# Calculating the Semi-Analytic FStat for detector H1
cmdline="$pfs_code $pfs_CL --DataFiles='./testSFTs/H-1_H1*'"
#echo $cmdline

if ! resPFSH1=`eval $cmdline 2> /dev/null`; then
    echo "Error ... something failed running '$pfs_code' ..."
    exit 1;
fi

res2PFSH1=`echo $resPFSH1 | awk '{printf "%.1f", $1}'`

#echo "$res2PFSH1"
############################################################
# Calculating the Semi-Analytic FStat for detector H2
cmdline="$pfs_code $pfs_CL --DataFiles='./testSFTs/H-1_H2*'"
#echo $cmdline

if ! resPFSH2=`eval $cmdline 2> /dev/null`; then
    echo "Error ... something failed running '$pfs_code' ..."
    exit 1;
fi

res2PFSH2=`echo $resPFSH2 | awk '{printf "%.1f", $1}'`

############################################################
# Calculating the Semi-Analytic FStat for detector L1
cmdline="$pfs_code $pfs_CL --DataFiles='./testSFTs/L-1_L*'"
#echo $cmdline

if ! resPFSL1=`eval $cmdline 2> /dev/null`; then
    echo "Error ... something failed running '$pfs_code' ..."
    exit 1;
fi

res2PFSL1=`echo $resPFSL1 | awk '{printf "%.1f", $1}'`

############################################################
# Calculating the Semi-Analytic FStat for detector G1
cmdline="$pfs_code $pfs_CL --DataFiles='./testSFTs/G-1_G*'"
#echo $cmdline

if ! resPFSG1=`eval $cmdline 2> /dev/null`; then
    echo "Error ... something failed running '$pfs_code' ..."
    exit 1;
fi

res2PFSG1=`echo $resPFSG1 | awk '{printf "%.1f", $1}'`

############################################################
# Calculating the Semi-Analytic FStat for detector H1+H2
cmdline="$pfs_code $pfs_CL --DataFiles='./testSFTs/H-1_H*'"
#echo $cmdline

if ! resPFSH1H2=`eval $cmdline 2> /dev/null`; then
    echo "Error ... something failed running '$pfs_code' ..."
    exit 1;
fi

res2PFSH1H2=`echo $resPFSH1H2 | awk '{printf "%.1f", $1}'`

############################################################
# Calculating the Semi-Analytic FStat for detector H1+H2+G1
cmdline="$pfs_code $pfs_CL --DataFiles='./testSFTs/'H-1_H*  G-1_G*''"
#echo $cmdline

if ! resPFSH1H2G1=`eval $cmdline 2> /dev/null`; then
    echo "Error ... something failed running '$pfs_code' ..."
    exit 1;
fi

res2PFSH1H2G1=`echo $resPFSH1H2G1 | awk '{printf "%.1f", $1}'`
##echo "res2PFSH1H2G1=$res2PFSH1H2G1"

############################################################
# Calculating the Semi-Analytic FStat for detector H1+H2+L1+G1
cmdline="$pfs_code $pfs_CL --DataFiles='./testSFTs/*-1*'"
#echo $cmdline

if ! resPFSH1H2L1G1=`eval $cmdline 2> /dev/null`; then
    echo "Error ... something failed running '$pfs_code' ..."
    exit 1;
fi

res2PFSH1H2L1G1=`echo $resPFSH1H2L1G1 | awk '{printf "%.1f", $1}'`

############################################################

#echo
#echo "-------------------------------------------------------------------------"
#echo " STEP 3: run CFS_v2 with perfect match, for single detector ($IFO1 and $IFO3)"
#echo "-------------------------------------------------------------------------"

## common cmdline-options for v2
cfs_CL="--Freq=$f0 --Alpha=$Alpha --Delta=$Delta --f1dot=$f1dot --TwoFthreshold=0 --refTime=$refTime"


############################################################
cmdline="$cfsv2_code $cfs_CL --DataFiles='$SFTdir/H-1_H1*'";
#echo $cmdline;
if ! resCFSH1=`eval $cmdline 2> /dev/null`; then
    echo "Error ... something failed running '$pfs_code' ..."
    exit 1;
fi

res2CFSH1=`echo $resCFSH1 | awk '{printf "%.1f", $1}'`
echo "$res2CFSH1"

############################################################
cmdline="$cfsv2_code $cfs_CL --DataFiles='$SFTdir/H-1_H2*'";
#echo $cmdline;
if ! resCFSH2=`eval $cmdline 2> /dev/null`; then
    echo "Error ... something failed running '$pfs_code' ..."
    exit 1;
fi

res2CFSH2=`echo $resCFSH2 | awk '{printf "%.1f", $1}'`
############################################################
cmdline="$cfsv2_code $cfs_CL --DataFiles='$SFTdir/L-1_L1*'";
#echo $cmdline;
if ! resCFSL1=`eval $cmdline 2> /dev/null`; then
    echo "Error ... something failed running '$pfs_code' ..."
    exit 1;
fi

res2CFSL1=`echo $resCFSL1 | awk '{printf "%.1f", $1}'`
############################################################

cmdline="$cfsv2_code $cfs_CL --DataFiles='$SFTdir/G-1_G1*'";
#echo $cmdline;
if ! resCFSG1=`eval $cmdline 2> /dev/null`; then
    echo "Error ... something failed running '$pfs_code' ..."
    exit 1;
fi

res2CFSG1=`echo $resCFSG1 | awk '{printf "%.1f", $1}'`
############################################################

cmdline="$cfsv2_code $cfs_CL --DataFiles='./testSFTs/H-1_H*'"
#echo $cmdline

if ! resCFSH1H2=`eval $cmdline 2> /dev/null`; then
    echo "Error ... something failed running '$pfs_code' ..."
    exit 1;
fi

res2CFSH1H2=`echo $resCFSH1H2 | awk '{printf "%.1f", $1}'`
############################################################

cmdline="$cfsv2_code $cfs_CL --DataFiles='./testSFTs/*-1*'"
#echo $cmdline

if ! resCFSH1H2L1G1=`eval $cmdline 2> /dev/null`; then
    echo "Error ... something failed running '$pfs_code' ..."
    exit 1;
fi

res2CFSH1H2L1G1=`echo $resCFSH1H2L1G1 | awk '{printf "%.1f", $1}'`
############################################################

echo " $Alpha    $Delta    $aPlus    $aCross    $psi     $phi0    $f0   $f1dot    $res2PFSH1   $res2CFSH1" >> resH1rand.txt
echo " $Alpha    $Delta    $aPlus    $aCross    $psi     $phi0    $f0   $f1dot    $res2PFSH2   $res2CFSH2" >> resH2rand.txt
echo " $Alpha    $Delta    $aPlus    $aCross    $psi     $phi0    $f0   $f1dot    $res2PFSL1   $res2CFSL1" >> resL1rand.txt
echo " $Alpha    $Delta    $aPlus    $aCross    $psi     $phi0    $f0   $f1dot    $res2PFSG1   $res2CFSG1" >> resG1rand.txt
echo " $Alpha    $Delta    $aPlus    $aCross    $psi     $phi0    $f0   $f1dot    $res2PFSH1H2   $res2CFSH1H2" >> resH1H2rand.txt
echo " $Alpha    $Delta    $aPlus    $aCross    $psi     $phi0    $f0   $f1dot    $res2PFSH1H2L1G1   $res2CFSH1H2L1G1" >> resH1H2L1G1rand.txt
    iteration=$(($iteration + 1));
done
## end of while-loop

