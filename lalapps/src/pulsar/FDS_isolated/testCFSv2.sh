#!/bin/bash

## run all LALApps programs with memory debugging
export LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},memdbg"

## make sure we work in 'C' locale here to avoid awk sillyness
LC_ALL_old=$LC_ALL
export LC_ALL=C

builddir="./";
injectdir="../Injections/"

## ----- user-controlled level of debug-output detail
if [ -n "$DEBUG" ]; then
    debug=${DEBUG}
else
    debug=0	## default=quiet
fi

##---------- names of codes and input/output files
mfd_code="${injectdir}lalapps_Makefakedata_v4"
saf_code="${builddir}lalapps_SemiAnalyticF"
cfs_code="${builddir}lalapps_ComputeFStatistic"
cmp_code="${builddir}lalapps_compareFstats"
## allow user to specify a different CFSv2 version to test by passing as cmdline-argument
if test $# -eq 0 ; then
    cfsv2_code="${builddir}lalapps_ComputeFStatistic_v2"
else
    cfsv2_code="$@"
fi

SFTdir="./testCFSv2_sfts"

if [ -n "${LALPULSAR_DATADIR}" ]; then
    mfd_code="${mfd_code} -E ${LALPULSAR_DATADIR}"
    saf_code="${saf_code} -E ${LALPULSAR_DATADIR}"
    cfs_code="${cfs_code} -E ${LALPULSAR_DATADIR}"
    cfsv2_code="${cfsv2_code} -E ${LALPULSAR_DATADIR}"
else
    echo
    echo "Need environment-variable LALPULSAR_DATADIR to be set to"
    echo "your ephemeris-directory (e.g. /usr/local/share/lalpulsar)"
    echo "This might indicate an incomplete LAL+LALPULSAR installation"
    echo
    exit 1
fi

## without noise-weights, CFSv1 and CFSv2 should agree extremely well,
Ftolerance_NWoff=0.01
## with noise-weights to calculations differ, so we need much more tolerance
Ftolerance_NWon=0.1

Dterms=8
# ---------- fixed parameter of our test-signal
Tsft=1800;
startTime=711595934
duration=144000		## 40 hours

mfd_FreqBand=2.0;

Alpha=2.0
Delta=-0.5

h0=1
cosi=-0.3

psi=0.6
phi0=1.5

Freq=100.12345
f1dot=-1e-10;

## mfd-specific bands
mfd_fmin=$(echo $Freq $mfd_FreqBand | awk '{printf "%g", $1 - $2 / 2.0}');

## cfs search bands
NFreq=500;
cfs_FreqBand=$(echo $duration | awk '{printf "%.16g", 1.0 / $1 }');	## fix band to 1/T so we're close to signal peak always
cfs_Freq=$(echo $Freq $cfs_FreqBand | awk '{printf "%.16g", $1 - $2 / 2.0}');
cfs_dFreq=$(echo $cfs_FreqBand $NFreq | awk '{printf "%.16g", $1 / $2 }');
cfs_nCands=$NFreq	## toplist length: keep all cands

## unfortunately CFSv1 has a different band-convention resulting in one more frequency-bin
## so we compensate for that by inputting a slightly smaller band in CFSv1 to get the same
## bins for comparisong
cfs_FreqBand_v1=$(echo $cfs_FreqBand $cfs_dFreq | awk '{printf "%g", $1 - 0.5 * $2}' )

cfs_f1dotBand=0;
cfs_f1dot=$(echo $f1dot $cfs_f1dotBand | awk '{printf "%.16g", $1 - $2 / 2.0}');
##Nf1dot=10
cfs_df1dot=1 ##$(echo $cfs_f1dotBand $Nf1dot | awk '{printf "%g", $1 / $2}');

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
mfd_CL="${saf_CL} --fmin=$mfd_fmin --Band=$mfd_FreqBand --Freq=$Freq --outSFTbname=$SFTdir/testSFT --f1dot=$f1dot --outSFTv1"
if [ "$haveNoise" = true ]; then
    mfd_CL="$mfd_CL --noiseSqrtSh=$sqrtSh";
fi

cmdline="$mfd_code $mfd_CL --randSeed=1"
echo $cmdline;
if ! eval "$cmdline &> /dev/null"; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

echo
echo -n "Running '$saf_code' ... "
cmdline="$saf_code $saf_CL --sqrtSh=$sqrtSh"
echo $cmdline
if ! resF=`eval "$cmdline  2> /dev/null"`; then
    echo "Error ... something failed running '$saf_code' ..."
    exit 1;
fi
res2F=$(echo $resF | awk '{printf "%g", 2.0 * $1}')
echo "The SemiAnalyticF calculations predicts: 2F = $res2F"

echo
echo "----------------------------------------------------------------------"
echo "STEP 2: run CFS_v1 with perfect match"
echo "----------------------------------------------------------------------"
echo
outfile_v1="Fstat_v1.dat";
## common cmdline-options for v1 and v2
cfs_CL="--IFO=$IFO --Alpha=$Alpha --Delta=$Delta --Freq=$cfs_Freq --dFreq=$cfs_dFreq --f1dot=$cfs_f1dot --f1dotBand=$cfs_f1dotBand --df1dot=$cfs_df1dot --DataFiles='$SFTdir/testSFT*' --Dterms=${Dterms} --NumCandidatesToKeep=${cfs_nCands}"
if [ "$haveNoise" = false ]; then
    cfs_CL="$cfs_CL --SignalOnly"
fi

cmdline="$cfs_code $cfs_CL  --outputFstat=$outfile_v1 --expLALDemod=0 --Fthreshold=0 --FreqBand=$cfs_FreqBand_v1"
echo $cmdline;

if ! eval "$cmdline &> /dev/null"; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo " STEP 3: run CFS_v2 with perfect match"
echo "----------------------------------------------------------------------"
echo
outfile_v2NWon="Fstat_v2NWon.dat";
cmdlineNoiseWeightsOn="$cfsv2_code $cfs_CL --outputFstat=$outfile_v2NWon --TwoFthreshold=0 --FreqBand=$cfs_FreqBand --UseNoiseWeights=true"
echo $cmdlineNoiseWeightsOn;
if ! eval "$cmdlineNoiseWeightsOn &> /dev/null"; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1;
fi

outfile_v2NWoff="Fstat_v2NWoff.dat";
cmdlineNoiseWeightsOff="$cfsv2_code $cfs_CL --outputFstat=$outfile_v2NWoff --TwoFthreshold=0 --FreqBand=$cfs_FreqBand --UseNoiseWeights=false"
echo $cmdlineNoiseWeightsOff;
if ! eval "$cmdlineNoiseWeightsOff &> /dev/null"; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1;
fi


echo
echo "----------------------------------------"
echo " STEP 4: Comparing results: "
echo "----------------------------------------"

## work around toplist-sorting bugs in CFSv2: manually sort before comparing
sort $outfile_v1 > __tmp_sorted && mv __tmp_sorted $outfile_v1
sort $outfile_v2NWoff > __tmp_sorted && mv __tmp_sorted $outfile_v2NWoff
sort $outfile_v2NWon > __tmp_sorted && mv __tmp_sorted $outfile_v2NWon

echo
cmdline="$cmp_code -1 ./$outfile_v1 -2 ./$outfile_v2NWoff --clusterFiles=0 --Ftolerance=$Ftolerance_NWoff"
echo -n $cmdline
if ! eval $cmdline; then
    echo "==> OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "	==> OK."
fi

cmdline="$cmp_code -1 ./$outfile_v1 -2 ./$outfile_v2NWon --clusterFiles=0 --Ftolerance=$Ftolerance_NWon"
echo -n $cmdline
if ! eval $cmdline; then
    echo "==> OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "	==> OK."
fi
echo

## -------------------------------------------
## clean up files
## -------------------------------------------
if [ -z "$NOCLEANUP" ]; then
    rm -rf $SFTdir $outfile_v1 $outfile_v2NWon $outfile_v2NWoff Fstats Fstats.log
fi

## restore original locale, just in case someone source'd this file
export LC_ALL=$LC_ALL_old
