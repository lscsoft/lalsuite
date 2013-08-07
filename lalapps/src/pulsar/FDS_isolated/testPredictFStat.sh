#!/bin/sh

## run all LALApps programs with memory debugging
export LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},memdbg"

## test PredictFStat by comparison with SemiAnalyticF
extra_args="$@"

builddir="./";
injectdir="../Injections/"

##---------- names of codes and input/output files
mfd_code="${injectdir}lalapps_Makefakedata_v4"
saf_code="${builddir}lalapps_SemiAnalyticF"
pfs_code="${builddir}lalapps_PredictFStat"

SFTdir="./testPredictFStat_sfts"

if [ -n "${LALPULSAR_DATADIR}" ]; then
    mfd_code="${mfd_code} -E ${LALPULSAR_DATADIR}"
    saf_code="${saf_code} -E ${LALPULSAR_DATADIR}"
    pfs_code="${pfs_code} -E ${LALPULSAR_DATADIR}"
else
    echo
    echo "Need environment-variable LALPULSAR_DATADIR to be set to"
    echo "your ephemeris-directory (e.g. /usr/local/share/lalpulsar)"
    echo "This might indicate an incomplete LAL+LALPULSAR installation"
    echo
    exit 1
fi

Ftolerance=0.05
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
mfd_fmin=$(echo $Freq $mfd_FreqBand | awk '{printf "%g", $1 - $2 / 2.0}');

noiseSqrtSh=1

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
mfd_CL="${saf_CL} --fmin=$mfd_fmin --Band=$mfd_FreqBand --Freq=$Freq --outSFTbname=$SFTdir --noiseSqrtSh=$noiseSqrtSh --randSeed=1"

cmdline="$mfd_code $mfd_CL";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

echo
cmdline="$saf_code $saf_CL  --sqrtSh=$noiseSqrtSh"
echo $cmdline
if ! tmp=`eval $cmdline 2> /dev/null`; then
    echo "Error ... something failed running '$saf_code' ..."
    exit 1;
fi
echo  "ok."
resSAF=`echo $tmp | awk '{printf "%g", 2.0 * $1}'`

echo
echo "----------------------------------------------------------------------"
echo "STEP 2: run PFS once with noise-SFTs"
echo "----------------------------------------------------------------------"
echo
outfile_pfs="__tmp_PFS.dat";

## common cmdline-options for v1 and v2
pfs_CL=" --Alpha=$Alpha --Delta=$Delta --IFO=$IFO --h0=$h0 --cosi=$cosi --psi=$psi --Freq=$Freq --DataFiles='$SFTdir/*' --outputFstat=$outfile_pfs"
cmdline="$pfs_code $pfs_CL"
echo $cmdline
if ! tmp=`eval $cmdline`; then
    echo "Error.. something failed when running '$pfs_code' ..."
    exit 1
fi
resPFS2=`echo $tmp | awk '{printf "%g", $1}'`

echo
echo "----------------------------------------------------------------------"
echo " STEP 3: run PFS without using noise-floor estimation (assume Sh=1)"
echo "----------------------------------------------------------------------"
echo
outfile_pfs2="__tmp_PFS2.dat";

## common cmdline-options for v1 and v2
pfs_CL=" --Alpha=$Alpha --Delta=$Delta --IFO=$IFO --h0=$h0 --cosi=$cosi --psi=$psi --Freq=$Freq --DataFiles='$SFTdir/*' --outputFstat=$outfile_pfs2 --SignalOnly"
cmdline="$pfs_code $pfs_CL"
echo $cmdline
if ! tmp=`eval $cmdline`; then
    echo "Error.. something failed when running '$pfs_code' ..."
    exit 1
fi
resPFS1=`echo $tmp | awk '{printf "%g", $1}'`

echo
echo "----------------------------------------"
echo " STEP 4: Comparing results: "
echo "----------------------------------------"
echo
echo "SemiAnalyticF:        2F_SA  = $resSAF"
echo "PredictFStat [Sh=1]:  2F_PF1 = $resPFS1"
echo "PredictFStat [SFTs]:  2F_PF2 = $resPFS2"
echo

eps1=$(echo $resSAF $resPFS1 | awk '{printf "%d", int ( 1000 * sqrt( ($1 - $2)*($1 - $2) ) / $1 )}');
eps2=$(echo $resSAF $resPFS2 | awk '{printf "%d", int ( 1000 * sqrt( ($1 - $2)*($1 - $2) ) / $1 )}');

res=0;
if test $eps1 -gt 10; then
    echo "==> relative deviation F_PF1 wrt F_SA: $eps1/1000 larger than 1% ... FAILED."
    echo
    res=1;
else
    echo "==> relative deviation F_PF1 wrt F_SA: $eps1/1000 smaller than 1% ... PASSED."
fi

if test $eps2 -gt 150; then
    echo "==> relative deviation F_PF2 wrt F_SA: $eps2/1000 larger than 15% ... FAILED."
    echo
    res=1;
else
    echo "==> relative deviation F_PF2 wrt F_SA: $eps2/1000 smaller than 15% ... PASSED."
    echo
fi

## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm -rf $SFTdir $outfile_pfs $outfile_pfs2;
fi

exit $res;

