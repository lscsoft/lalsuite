#!/bin/sh
## A simple test script to compare the F-Statistic of a single and two detectors.
## This uses ComputeFStatistic_v2 to calcultae the F-Statistic and the makefakedata 
## to generate v2 version of SFT.

## Author: Iraj Gholami

## take user-arguments for CFS-v2:
extra_args="$@"

##---------- names of codes and input/output files  
saf_code="lalapps_SemiAnalyticF"
mfd_code="lalapps_Makefakedata"
cfs_code="lalapps_ComputeFStatistic"
cfsv2_code="ComputeFStatistic_v2"
cmp_code="compareFstats"
pfs_code="lalapps_PredictFStat"

SFTdir="./testSFTs"

# test if LAL_DATA_PATH has been set ... needed to locate ephemeris-files
if [ -z "$LAL_DATA_PATH" ]; then
    echo
    echo "Need environment-variable LAL_DATA_PATH to point "
    echo "to your ephemeris-directory (e.g. /usr/local/share/lal)"
    if [ -n "$LAL_PREFIX" ]; then
	echo "You have LAL_PREFIX set, I suggest setting 'LAL_DATA_PATH=\$LAL_PREFIX/share/lal'"
    fi
    echo
    exit 1
fi

# ---------- fixed parameter of our test-signal
##Tsft=60;
Tsft=1800;
startTime=711595933
refTime=711595933     ##$startTime 
duration=144000	      ## 27.7777 hours

mfd_FreqBand=2.0;

Alpha=2.0
Delta=-0.5

aPlus=1.0e-20
aCross=0.4e-20
##aPlus=0
##aCross=0

psi=0
phi0=0

freq=100.0
mfd_fmin=$(echo $freq $mfd_FreqBand | awk '{printf "%g", $1 - $2 / 2.0}');

echo "mfd_fmin = $mfd_fmin"

f1dot=1e-8
df1dot=0.3e-8	## search about 3 spindown-values

##noiseSqrtSh=3e-23
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
saf_CL="--latitude=$Delta  --longitude=$Alpha --detector=$IFO1 --Tsft=$Tsft --startTime=$startTime --duration=$duration --aPlus=$aPlus --aCross=$aCross --psi=$psi --phi0=$phi0"
# concatenate this with the mfd-specific switches:
mfd_CL="${saf_CL} --fmin=$mfd_fmin --Band=$mfd_FreqBand --f0=$freq --outSFTbname=$SFTdir/ --f1dot=$f1dot  --refTime=$refTime --noiseSqrtSh=$sqrtSh"
## --outSFTv1"

# this part of the command-line is compatible with SemiAnalyticF:
saf2_CL="--latitude=$Delta  --longitude=$Alpha --detector=$IFO2 --Tsft=$Tsft --startTime=$startTime --duration=$duration --aPlus=$aPlus --aCross=$aCross --psi=$psi --phi0=$phi0"
# concatenate this with the mfd-specific switches:
mfd2_CL="${saf2_CL} --fmin=$mfd_fmin --Band=$mfd_FreqBand --f0=$freq --outSFTbname=$SFTdir/ --f1dot=$f1dot  --refTime=$refTime --noiseSqrtSh=$sqrtSh"
## --outSFTv1"

cmdline="$mfd_code $mfd_CL";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

echo
cmdline="$mfd_code $mfd2_CL";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo " STEP 2: Calculate F-Statistic Semi-Analytically"
echo "----------------------------------------------------------------------"
echo 

cmdline="$saf_code $saf_CL --sqrtSh=$sqrtSh"	
echo $cmdline
if ! resF=`eval $cmdline 2> /dev/null`; then
    echo "Error ... something failed running '$saf_code' ..."
    exit 1;
fi
#echo  "ok."
res2F=`echo $resF | awk '{printf "%g", 2.0 * $1}'`
echo

# this part of the command-line is compatible with PredictFStat:
pfs_CL="--aPlus=$aPlus --aCross=$aCross --psi=$psi --phi=$phi0 --Freq=100.0 --DataFiles='./testSFTs/H-1_H1*' --Delta=$Delta --Alpha=$Alpha"

cmdline="$pfs_code $pfs_CL"
echo $cmdline

if ! resPFS1=`eval $cmdline 2> /dev/null`; then
    echo "Error ... something failed running '$pfs_code' ..."
    exit 1;
fi
#echo  "ok."

res2PFS1=`echo $resPFS1 | awk '{printf "%g", $1}'`
echo

pfs2_CL="--aPlus=$aPlus --aCross=$aCross --psi=$psi --phi=$phi0 --Freq=100.0 --DataFiles='./testSFTs/H-1_H*' --Delta=$Delta --Alpha=$Alpha"

cmdline="$pfs_code $pfs2_CL"
echo $cmdline

if ! resPFS2=`eval $cmdline 2> /dev/null`; then
    echo "Error ... something failed running '$pfs_code' ..."
    exit 1;
fi
#echo  "ok."

res2PFS2=`echo $resPFS2 | awk '{printf "%g", $1}'`
echo

echo "The SemiAnalyticF for $IFO1,    predicts: 2F = $res2F"
echo "The PredictFStat  for $IFO1,    predicts: 2F = $res2PFS1"
echo "The PredictFStat  for $IFO1+$IFO2, predicts: 2F = $res2PFS2"

echo 
echo "----------------------------------------------------------------------"
echo " STEP 3: run CFS_v2 with perfect match, for single detector"
echo "----------------------------------------------------------------------"
echo

## cmdline-options for v2, single detector   
cfs1_CL="--Freq=$freq --Alpha=$Alpha --Delta=$Delta --f1dot=$f1dot --f1dotBand=$f1dot --df1dot=$df1dot --Fthreshold=0 --DataFiles='$SFTdir/H-1_H1*' --refTime=$refTime $whatNoise"
    
cmdline="$cfsv2_code $cfs1_CL --outputFstat=Fstat_v2_1.dat";
echo $cmdline;

if ! eval time $cmdline; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1
fi

echo    
echo "----------------------------------------------------------------------"
echo " STEP 4: run CFS_v2 with perfect match, for two detectors"
echo "----------------------------------------------------------------------"
echo

## cmdline-options for v2, two detectors  
cfs2_CL="--Freq=$freq --Alpha=$Alpha --Delta=$Delta --f1dot=$f1dot --f1dotBand=$f1dot --df1dot=$df1dot --Fthreshold=0 --DataFiles='$SFTdir/H-1_H*' --refTime=$refTime $whatNoise"

cmdline="$cfsv2_code $cfs2_CL --outputFstat=Fstat_v2_2.dat $extra_args";
echo $cmdline;

if ! eval time $cmdline; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1;
fi

echo
echo "----------------------------------------"
echo " STEP 5: Comparing results: "
echo "----------------------------------------"
echo
echo "Fstat_v2_1.dat: "
cat Fstat_v2_1.dat
echo "Fstat_v2_2.dat: "
cat Fstat_v2_2.dat

echo
##cmdline="$cmp_code -1 ./Fstat_v2_1.dat -2 ./Fstat_v2_2.dat --clusterFiles=0 --Ftolerance=0.1"
##echo $cmdline
##if ! eval $cmdline; then
##    echo "OUCH... files differ. Something might be wrong..."
##    exit 2
##else
##    echo "OK."
##fi
##
##echo
