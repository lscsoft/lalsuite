#!/bin/sh

## take extra user-arguments
extra_args="$@"

## allow 'make test' to work from builddir != srcdir
if [ -n "${srcdir}" ]; then
    builddir="./";
    injectdir="../Injections/"
else
    srcdir=./
    builddir="~/Software/development/lalapps-BUILD-dbg/src/pulsar/FDS_isolated/"
fi

##---------- names of codes and input/output files
mfd_code="${injectdir}lalapps_Makefakedata_v4"
lft_code="${builddir}lalapps_computeLFTfromSFTs"
dumpSFT_code="${injectdir}lalapps_dumpSFT"
extract_code="lalapps_ConvertToSFTv2"

SFTdir1="./testSFTs"
SFTdir2="./testSFTs2"

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

# ---------- fixed parameter of our test-signal
Tsft=1800;
startTime=711595934
refTime=701595833  ## $startTime
duration=144000		## 40 hours

Alpha=2.0
Delta=-0.5

h0=1
cosi=-0.3

psi=0.6
phi0=1.5

Freq=101.12345

mfd_FreqBand=2
mfd_fmin=100.1234
mfd_fmax=102.1234

f1dot=-1e-10;

noiseSqrtSh=0

TDDfile="tdd.dat"
LFTfile="out.lft"
LFTascii="lft.dat"

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
if [ ! -d "$SFTdir1" ]; then
    mkdir $SFTdir1;
else
    rm -f $SFTdir1/*;
fi
if [ ! -d "$SFTdir2" ]; then
    mkdir $SFTdir2;
else
    rm -f $SFTdir2/*;
fi

# this part of the command-line is compatible with SemiAnalyticF:
mfd_CL1=" --Alpha=$Alpha --Delta=$Delta --IFO=$IFO --Tsft=$Tsft --h0=$h0 --cosi=$cosi --psi=$psi --phi0=$phi0 --fmin=$mfd_fmin --Band=$mfd_FreqBand --Freq=$Freq --f1dot=$f1dot --refTime=$refTime --TDDfile=$TDDfile --timestampsFile=./ts.dat --generationMode=1 --outSFTbname=$SFTdir1/"

mfd_CL2=" --Alpha=$Alpha --Delta=$Delta --IFO=$IFO --Tsft=$Tsft --h0=$h0 --cosi=$cosi --psi=$psi --phi0=$phi0 --fmin=0 --Band=103 --Freq=$Freq --f1dot=$f1dot --refTime=$refTime --timestampsFile=./ts.dat --generationMode=0 --outSFTbname=${SFTdir2}/"

if [ "$haveNoise" = true ]; then
    mfd_CL1="$mfd_CL1 --noiseSqrtSh=$sqrtSh";
    mfd_CL2="$mfd_CL2 --noiseSqrtSh=$sqrtSh";
fi

cmdline="$mfd_code $mfd_CL1";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

## ---------- make high-freq SFT then extract sub-band
cmdline="$mfd_code $mfd_CL2";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

## extract relevant frequency-band
extract_CL="--inputSFTs='$SFTdir2/*.sft' --fmin=$mfd_fmin --fmax=$mfd_fmax --outputDir=$SFTdir2 --descriptionMisc=Band"
cmdline="$extract_code $extract_CL"
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$extract_code' ..."
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo "STEP 2: run computeLFTfromSFTs"
echo "----------------------------------------------------------------------"
echo
lft_CL=" --inputSFTs='$SFTdir1/*.sft' --outputLFT=$LFTfile --upsampling=1 2> tddSFT.dat"  ##  --fmin=100.0123451 --fmax=100.1994

cmdline="$lft_code $lft_CL";
echo $cmdline;

if ! eval $cmdline; then
    echo "Error.. something failed when running '$lft_code' ..."
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo " STEP 3: convert output LFT into ASCII"
echo "----------------------------------------------------------------------"
echo
cmdline="$dumpSFT_code --SFTfiles=$LFTfile --noHeader > $LFTascii";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$dumpSFT_code' ..."
    exit 1;
fi

## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm -rf $SFTdir1 $SFTdir2
fi
