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
cfs_code="${builddir}lalapps_CFSv3"
dumpSFT_code="${injectdir}lalapps_dumpSFT"
extract_code="lalapps_ConvertToSFTv2"

noiseSFTdir="./noiseSFTs"
SFTdir="./SFTdir"


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
refTime=701595833  ## $startTime

Alpha=2.0
Delta=-0.5

h0=1
cosi=-0.3

psi=0.6
phi0=1.5


fmin=9.12345
Freq=12.12345
fmax=15.12277    ## make sure we get an effective band of exactly 6Hz
fUpper=20

f1dot=-1e-10;

noiseSqrtSh=0

TDDfile="tdd.dat"
LFTfile="out.lft"
LFTascii="lft.dat"

## ------------------------------------------------------------
IFO=H1

## Prepare SFT directories
if [ ! -d "$SFTdir" ]; then
    mkdir $SFTdir;
else
    rm -f $SFTdir/*;
fi


##--------------------------------------------------
## test starts here
##--------------------------------------------------

if [ "$noiseSqrtSh" != 0 ]; then
    echo
    echo "----------------------------------------------------------------------"
    echo " STEP 0: Generate (full-band) noise SFTs"
    echo "----------------------------------------------------------------------"
    echo

    if [ ! -d "$noiseSFTdir" ]; then
	mkdir $noiseSFTdir;
    else
	rm -f $noiseSFTdir/*;
    fi

    mfd_CL=" --IFO=$IFO --Tsft=$Tsft --h0=0 --fmin=0 --Band=$fUpper --timestampsFile=./ts.dat --generationMode=0 --outSFTbname=${noiseSFTdir}/ --noiseSqrtSh=$noiseSqrtSh"

    cmdline="$mfd_code $mfd_CL";
    echo $cmdline;
    if ! eval $cmdline; then
	echo "Error.. something failed when running '$mfd_code' ..."
	exit 1
    fi

fi ## if noiseSqrtSh


echo
echo "----------------------------------------------------------------------"
echo " STEP 1: Generate (full-band) timeseries and SFTs"
echo "----------------------------------------------------------------------"
echo

mfd_CL=" --lineFeature --Alpha=$Alpha --Delta=$Delta --IFO=$IFO --Tsft=$Tsft --h0=$h0 --cosi=$cosi --psi=$psi --phi0=$phi0 --fmin=0 --Band=$fUpper --Freq=$Freq --f1dot=$f1dot --refTime=$refTime --TDDfile=$TDDfile --timestampsFile=./ts.dat --generationMode=0 --outSFTbname=${SFTdir}/ -v1"

if [ "$noiseSqrtSh" != 0 ]; then
    mfd_CL="$mfd_CL --noiseSFTs='${noiseSFTdir}/*.sft'"
fi

cmdline="$mfd_code $mfd_CL";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo "STEP 2: run CFSv3 on small band extracted from SFTs"
echo "----------------------------------------------------------------------"
echo
cfs_CL=" --inputSFTs='${SFTdir}/*.sft' --outputLFT=$LFTfile --fmin=$fmin --fmax=$fmax --upsampling=1 2> tddSFT.dat"

cmdline="$cfs_code $cfs_CL";
echo $cmdline;

if ! eval $cmdline; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1
fi

## convert output LFT into ASCII"
cmdline="$dumpSFT_code --SFTfiles=$LFTfile --noHeader > $LFTascii";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$dumpSFT_code' ..."
    exit 1;
fi


## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm -rf $SFTdir $noiseSFTdir
fi
