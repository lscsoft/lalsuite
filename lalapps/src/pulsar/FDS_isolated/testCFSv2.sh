#!/bin/sh

## take user-arguments for CFS-v2:
extra_args="$@"

##---------- names of codes and input/output files  
saf_code="lalapps_SemiAnalyticF"
mfd_code="lalapps_Makefakedata"
cfs_code="lalapps_ComputeFStatistic"
cfsv2_code="ComputeFStatistic_v2"
cmp_code="compareFstats"

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
Tsft=1800;
startTime=711595934
refTime=701595833  ## $startTime
duration=144000		## 40 hours

mfd_FreqBand=2.0;

Alpha=2.0
Delta=-0.5

h0=1
cosi=-0.3

psi=0.6
phi0=1.5

Freq=100.12345
mfd_fmin=$(echo $Freq $mfd_FreqBand | awk '{printf "%g", $1 - $2 / 2.0}');

echo "mfd_fmin = $mfd_fmin"

f1dot=0
df1dot=0.3e-8	## search about 3 spindown-values

noiseSqrtSh=5

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
mfd_CL="${saf_CL} --fmin=$mfd_fmin --Band=$mfd_FreqBand --Freq=$Freq --outSFTbname=$SFTdir/testSFT --f1dot=$f1dot  --refTime=$refTime --outSFTv1"
if [ "$haveNoise" = true ]; then
    mfd_CL="$mfd_CL --noiseSqrtSh=$sqrtSh";
fi

cmdline="$mfd_code $mfd_CL";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

echo 
echo -n "Running '$saf_code' ... "
cmdline="$saf_code $saf_CL --sqrtSh=$sqrtSh"	
echo $cmdline
if ! resF=`eval $cmdline 2> /dev/null`; then
    echo "Error ... something failed running '$saf_code' ..."
    exit 1;
fi
echo  "ok."
res2F=`echo $resF | awk '{printf "%g", 2.0 * $1}'`
echo "The SemiAnalyticF calculations predicts: 2F = $res2F"
  
echo 
echo "----------------------------------------------------------------------"
echo "STEP 2: run CFS_v1 with perfect match"
echo "----------------------------------------------------------------------"
echo

## common cmdline-options for v1 and v2    
cfs_CL="--IFO=$IFO --Freq=$Freq --Alpha=$Alpha --Delta=$Delta --f1dot=$f1dot --f1dotBand=$f1dot --df1dot=$df1dot --DataFiles='$SFTdir/testSFT*' --refTime=$refTime"
if [ "$haveNoise" = false ]; then
    cfs_CL="$cfs_CL --SignalOnly"
fi
    
cmdline="$cfs_code $cfs_CL  --outputFstat=Fstat_v1.dat --expLALDemod=1 --Fthreshold=0";
echo $cmdline;

if ! eval time $cmdline; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1
fi

echo    
echo "----------------------------------------------------------------------"
echo " STEP 3: run CFS_v2 with perfect match"
echo "----------------------------------------------------------------------"
echo
cmdline="$cfsv2_code $cfs_CL --outputFstat=Fstat_v2.dat --TwoFthreshold=0 $extra_args";
cmdline="$cmdline --outputLoudest=loudest.dat";
echo $cmdline;

if ! eval time $cmdline; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1;
fi

echo
echo "----------------------------------------"
echo " STEP 4: Comparing results: "
echo "----------------------------------------"
echo
echo "Fstat_v1.dat: "
cat Fstat_v1.dat
echo "Fstat_v2.dat: "
cat Fstat_v2.dat

echo
cmdline="$cmp_code -1 ./Fstat_v1.dat -2 ./Fstat_v2.dat --clusterFiles=0 --Ftolerance=0.1"
echo $cmdline
if ! eval $cmdline; then
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi

echo
