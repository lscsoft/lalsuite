#!/bin/sh

## test PredictFStat by comparison with SemiAnalyticF
extra_args="$@"

## allow 'make test' to work from builddir != srcdir
if [ -n "${srcdir}" ]; then
    builddir="./";
    injectdir="../Injections/"
else
    srcdir=.
fi

##---------- names of codes and input/output files
mfd_code="${injectdir}lalapps_Makefakedata_v4"
saf_code="${builddir}lalapps_SemiAnalyticF"
pfs_code="${builddir}lalapps_PredictFStat"

SFTdir="./testSFTs"

# test if LAL_DATA_PATH has been set ... needed to locate ephemeris-files
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
echo -n "Running '$saf_code' ... "
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
echo -n "Running '$pfs_code' ... "
echo $cmdline
if ! tmp=`eval $cmdline`; then
    echo "Error.. something failed when running '$pfs_code' ..."
    exit 1
fi
resPFS=`echo $tmp | awk '{printf "%g", $1}'`

echo
echo "----------------------------------------------------------------------"
echo " STEP 3: run PFS without using noise-SFT data-files"
echo "----------------------------------------------------------------------"
echo



## clean up files
rm -rf $SFTdir $outfile_pfs

echo
echo "----------------------------------------"
echo " STEP 4: Comparing results: "
echo "----------------------------------------"
echo
echo "SemiAnalyticF:        2F = $resSAF"
echo "PredictFStat [SFTs]:  2F = $resPFS"
echo

eps=$(echo $resSAF $resPFS | awk '{printf "%d", int ( 1000 * sqrt( ($1 - $2)*($1 - $2) ) / $1 )}');

if test $eps -gt 150; then
    echo "==> relative deviation $eps/1000 larger than 15% ... FAILED."
    echo
    exit 1;
else
    echo "==> relative deviation $eps/1000 smaller than 15% ... PASSED."
    echo
    exit 0;
fi


