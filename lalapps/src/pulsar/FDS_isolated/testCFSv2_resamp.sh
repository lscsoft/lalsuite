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
cmp_code="${builddir}lalapps_compareFstats"

## allow user to specify a different CFSv2 version to test by passing as cmdline-argument
if test $# -eq 0 ; then
    cfs_code="${builddir}lalapps_ComputeFStatistic_v2"
else
    cfs_code="$@"
fi

SFTdir="./testCFSv2_resamp_sfts"

if [ -n "${LALPULSAR_DATADIR}" ]; then
    mfd_code="${mfd_code} -E ${LALPULSAR_DATADIR}"
    cfs_code="${cfs_code} -E ${LALPULSAR_DATADIR}"
else
    echo
    echo "Need environment-variable LALPULSAR_DATADIR to be set to"
    echo "your ephemeris-directory (e.g. /usr/local/share/lalpulsar)"
    echo "This might indicate an incomplete LAL+LALPULSAR installation"
    echo
    exit 1
fi

# ---------- fixed parameter of our test-signal
Tsft=1800;
startTime=711595934
duration=144000		## 40 hours
refTime=611595934	## ~3.5 years prior to startTime

mfd_FreqBand=2.0;

Alpha=2.0
Delta=-0.5
AlphaBand=0.1
DeltaBand=0.1
dAlpha=0.05
dDelta=0.05

h0=1
cosi=-0.3
psi=0.6
phi0=1.5

Freq=100.12345
f1dot=-1e-10;

## mfd-specific bands
mfd_fmin=$(echo $Freq $mfd_FreqBand | awk '{printf "%g", $1 - $2 / 2.0}');

## cfs search bands
Nfreq=500
cfs_dFreq=$(echo $duration | awk '{printf "%.16g", 1.0 / ( 2.0 * $1 ) }');		## 1/(2T) frequency resolution
cfs_FreqBand=$(echo $Nfreq $cfs_dFreq | awk '{printf "%.16g", $1 * $2 - 0.5*$2 }');	## band corresponding to fixed number of frequency bins
cfs_Freq=$(echo $Freq $cfs_FreqBand | awk '{printf "%.16g", $1 - $2 / 2.0}');		## center search band on signal frequency

Nf1dot=10
cfs_df1dot=$(echo $duration | awk '{printf "%g", 0.005 / ($1 * $1) }');			## 1/(T^2) resolution in f1dot
cfs_f1dotBand=$(echo $Nf1dot $cfs_df1dot | awk '{printf "%.16g", $1 * $2 - 0.5*$2 }');	## band corresponding to fixed number of f1dot bins
cfs_f1dot=$(echo $f1dot $cfs_f1dotBand | awk '{printf "%.16g", $1 - $2 / 2.0}');	## center f1dot band on signal f1dot

cfs_nCands=100000	## toplist length: keep N cands

noiseSqrtSh=5
## ------------------------------------------------------------
if [ "$noiseSqrtSh" != 0 ]; then
    haveNoise=true;
else
    haveNoise=false;
fi

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

mfd_CL1="--refTime=${refTime} --Alpha=$Alpha --Delta=$Delta --Tsft=$Tsft --startTime=$startTime --duration=$duration --h0=$h0 --cosi=$cosi --psi=$psi --phi0=$phi0"
mfd_CL="${mfd_CL1} --fmin=$mfd_fmin --Band=$mfd_FreqBand --Freq=$Freq --outSFTbname=$SFTdir --f1dot=$f1dot"
if [ "$haveNoise" = true ]; then
    mfd_CL="$mfd_CL --noiseSqrtSh=$noiseSqrtSh";
fi

cmdline="$mfd_code $mfd_CL --randSeed=1 --IFO=H1"
echo $cmdline;
if ! eval "$cmdline &> /dev/null"; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

cmdline="$mfd_code $mfd_CL --randSeed=2 --IFO=L1"
echo $cmdline;
if ! eval "$cmdline &> /dev/null"; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo "STEP 2: run directed CFS_v2 with LALDemod method"
echo "----------------------------------------------------------------------"
echo
cfs_CL=" --refTime=${refTime} --Alpha=$Alpha --Delta=$Delta  --AlphaBand=$AlphaBand --DeltaBand=$DeltaBand --dAlpha=$dAlpha --dDelta=$dDelta --Freq=$cfs_Freq --FreqBand=$cfs_FreqBand --dFreq=$cfs_dFreq --f1dot=$cfs_f1dot --f1dotBand=$cfs_f1dotBand --df1dot=${cfs_df1dot} --DataFiles='$SFTdir/*.sft' --NumCandidatesToKeep=${cfs_nCands}"
if [ "$haveNoise" = false ]; then
    cfs_CL="$cfs_CL --SignalOnly"
fi

outfile_LD="Fstat_LD.dat";
cmdline="$cfs_code $cfs_CL  --outputFstat=$outfile_LD --useResamp=0"
echo $cmdline;

if ! eval "$cmdline &> /dev/null"; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo " STEP 3: run directed CFS_v2 with resampling"
echo "----------------------------------------------------------------------"
echo
outfile_RS="Fstat_RS.dat";

cmdline="$cfs_code $cfs_CL  --outputFstat=$outfile_RS --useResamp=1"
echo $cmdline;
if ! eval "$cmdline &> /dev/null"; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1;
fi

echo "----------------------------------------"
echo " STEP 4: Comparing results: "
echo "----------------------------------------"

sort $outfile_LD > __tmp_sorted && mv __tmp_sorted $outfile_LD
sort $outfile_RS > __tmp_sorted && mv __tmp_sorted $outfile_RS

## compare absolute differences instead of relative, allow deviations of up to sigma=sqrt(8)~2.8
echo
cmdline="$cmp_code -1 ./${outfile_LD} -2 ./${outfile_RS} --clusterFiles=0 --sigFtolerance --Ftolerance=0.5"
echo -n $cmdline
if ! eval $cmdline; then
    echo "==> OUCH... files differ. Something might be wrong..."
    exit 2
else
    echo "	==> OK."
fi

## -------------------------------------------
## clean up files
## -------------------------------------------
if [ -z "$NOCLEANUP" ]; then
    rm -rf $SFTdir $outfile_LD $outfile_RS
fi

## restore original locale, just in case someone source'd this file
export LC_ALL=$LC_ALL_old
