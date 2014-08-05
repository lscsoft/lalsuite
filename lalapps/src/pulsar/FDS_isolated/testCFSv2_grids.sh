#!/bin/sh

## set LAL debug level
echo "Setting LAL_DEBUG_LEVEL=${LAL_DEBUG_LEVEL:-msglvl1,memdbg}"
export LAL_DEBUG_LEVEL

## allow 'make test' to work from builddir != srcdir
if [ -z "${srcdir}" ]; then
    srcdir=`dirname $0`
fi

## take user-arguments for CFS-v2:
extra_args="$@"

builddir="./";
injectdir="../Injections/"

## ----- user-controlled level of debug-output detail
if [ -n "$DEBUG" ]; then
    debug=${DEBUG}
else
    debug=0	## default=quiet
fi

##---------- names of codes and input/output files
mfdv4_code="${injectdir}lalapps_Makefakedata_v4"
mfdv5_code="${injectdir}lalapps_Makefakedata_v5"
cfsv2_code="${builddir}lalapps_ComputeFStatistic_v2"
cmp_code="${builddir}lalapps_compareFstats"

# ---------- fixed parameter of our test-signal
Tsft=1800;
startTime=711595934
duration=144000		## 40 hours

Alpha=2.0
AlphaBand=0.1
dAlpha=0.04
Delta=-0.5
DeltaBand=0.1
dDelta=0.04

h0=1
cosi=-0.3

psi=0.6
phi0=1.5

Freq=100.12345
FreqBand=0.1
dFreq=0.04
f1dot=-1e-13
f1dotBand=1e-14
df1dot=4e-15

mfd_FreqBand=4.0;
mfd_fmin=$(echo $Freq $mfd_FreqBand | awk '{printf "%g", $1 - $2 / 2.0}');

echo "mfd_fmin = $mfd_fmin"

noiseSqrtSh=0

## ------------------------------------------------------------

if [ "$noiseSqrtSh" != 0 ]; then
    sqrtSh=$noiseSqrtSh
    haveNoise=true;
else
    sqrtSh=1;	## for SemiAnalyticF signal-only case
    haveNoise=false;
fi

IFO=LLO

##--------------------------------------------------
## test starts here
##--------------------------------------------------

echo
echo "----------------------------------------------------------------------"
echo " STEP 1: Generate Fake SFTs"
echo "----------------------------------------------------------------------"
echo

## --- for grid types 0,1,2,3,6, generate 40 hours of contiguous SFTs

## create SFT directory
SFTdir_40h="./testCFSv2_grids_SFTs_40h"
rm -rf $SFTdir_40h
mkdir $SFTdir_40h

# build MFD_v4 command line
mfd_CL=" --Alpha=$Alpha --Delta=$Delta --IFO=$IFO --Tsft=$Tsft --startTime=$startTime --duration=$duration --h0=$h0 --cosi=$cosi --psi=$psi --phi0=$phi0"
mfd_CL="${mfd_CL} --fmin=$mfd_fmin --Band=$mfd_FreqBand --Freq=$Freq --outSFTbname=$SFTdir_40h/testSFT --f1dot=$f1dot --outSFTv1"
if [ "$haveNoise" = true ]; then
    mfd_CL="$mfd_CL --noiseSqrtSh=$sqrtSh";
fi

## generate SFTs
cmdline="$mfdv4_code $mfd_CL";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdv4_code' ..."
    exit 1
fi

## --- for grid types 8,9, generate 2 SFTs spaced 5 days apart

## create SFT directory
SFTdir_5d="./testCFSv2_grids_SFTs_5d"
rm -rf $SFTdir_5d
mkdir $SFTdir_5d

## create timestamps file
cat <<EOF >"${SFTdir_5d}/timestamps.txt"
800000000 0
800432000 0
EOF

## build MFD_v5 command line
mfd_CL="--outSingleSFT=true --outSFTdir=${SFTdir_5d} --IFOs=H1 --sqrtSX=1.0 --timestampsFiles=${SFTdir_5d}/timestamps.txt --fmin=100 --Band=1.0 --randSeed=12345"

## generate SFTs
cmdline="$mfdv5_code $mfd_CL";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error: something failed when running '$mfdv5_code' ..."
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo " STEP 2: run CFS_v2 for various grid-types"
echo "----------------------------------------------------------------------"
echo

## common arguments for grid types 0,1,2,3,6
sky_CL="--Alpha=$Alpha --AlphaBand=$AlphaBand --dAlpha=$dAlpha --Delta=$Delta --DeltaBand=$DeltaBand --dDelta=$dDelta"
spin_CL="--Freq=$Freq --FreqBand=$FreqBand --dFreq=$dFreq --f1dot=$f1dot --f1dotBand=$f1dotBand --df1dot=$df1dot"
cfs_CL="--IFO=$IFO --DataFiles='${SFTdir_40h}/testSFT*' --TwoFthreshold=0 --Dterms=16 --FstatMethod=DemodOptC $extra_args"
if [ "$haveNoise" = false ]; then
    cfs_CL="$cfs_CL --SignalOnly"
fi

## ----- grid=0 : flat grid
echo "CFSv2 using gridType=0:"
cmdline="$cfsv2_code $sky_CL $spin_CL $cfs_CL --gridType=0 --outputFstat=./testCFSv2_grid0.dat";
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cmdline' ..."
    exit 1
fi

## ----- grid=1 : isotropic
echo "CFSv2 using gridType=1:"
cmdline="$cfsv2_code $sky_CL $spin_CL $cfs_CL --gridType=1 --outputFstat=./testCFSv2_grid1.dat";
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cmdline' ..."
    exit 1
fi

## ----- grid=2 : metric grid
echo "CFSv2 using gridType=2:"
cmdline="$cfsv2_code $sky_CL $spin_CL $cfs_CL --gridType=2 --metricType=1 --metricMismatch=0.1 --outputFstat=./testCFSv2_grid2.dat";
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cmdline' ..."
    exit 1
fi

## ----- grid=3 : same metric grid, specified as 'grid-file' passed on the commandline
echo "CFSv2 using gridType=3:"
cmdline="$cfsv2_code $sky_CL $spin_CL $cfs_CL --gridType=3 --outputFstat=./testCFSv2_grid3.dat --gridFile='{2.000000953674316 -0.4704365134239197; 2.068386077880859 -0.4704365134239197; 2.099475860595703 -0.492608368396759; 2.000000953674316 -0.4204375147819519; 2.082211971282959 -0.4204375147819519}'";
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cmdline' ..."
    exit 1
fi

## ----- grid=6 : recompute metric grid with --gridFile option
echo "recompute gridType=2 with CFSv2 using gridType=6:"
## extract a 6-column gridFile from the previous result output-file
gridFile="testCFSv2_grid2_grid.dat";
awk_extract6='{printf "%s %s %s %s %s %s\n", $1, $2, $3, $4, $5, $6 }'
sed '/^%.*/d' ./testCFSv2_grid2.dat | awk "$awk_extract6" >${gridFile}

cmdline="$cfsv2_code $cfs_CL --gridType=6 --gridFile=./${gridFile} --outputFstat=./testCFSv2_grid6.dat";
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cmdline' ..."
    exit 1
fi

## --- grid=8 : lattice tiling grid, square parameter space
echo "CFSv2 using gridType=8:"
cmdline="$cfsv2_code --Alpha=6.1 --Delta=1.2 --Freq=100.4 --FreqBand=5e-4 --f1dot=-1e-10 --f1dotBand=1e-10 --DataFiles='${SFTdir_5d}/*.sft' --TwoFthreshold=0 --gridType=8 --metricMismatch=0.5 --outputFstat=./testCFSv2_grid8.dat"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cmdline' ..."
    exit 1
fi

## --- grid=9 : lattice tiling grid, age-spindown-index parameter space
echo "CFSv2 using gridType=9:"
cmdline="$cfsv2_code --Alpha=6.1 --Delta=1.2 --Freq=100.4 --FreqBand=8e-5 --spindownAge=1e11 --minBraking=2 --maxBraking=5 --DataFiles='${SFTdir_5d}/*.sft' --TwoFthreshold=0 --gridType=9 --metricMismatch=0.5 --outputFstat=./testCFSv2_grid9.dat"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cmdline' ..."
    exit 1
fi

echo
echo "----------------------------------------"
echo " STEP 3: Compare to reference results: "
echo "----------------------------------------"
echo

for n in 0 1 2 3 6 8 9; do

    ## compare results
    echo "Comparing gridType=${n}:"
    cmdline="$cmp_code -1 ./testCFSv2_grid${n}.dat -2 ${srcdir}/testCFSv2_grid${n}.dat.ref.gz";
    echo $cmdline
    if ! eval $cmdline; then
        echo "OUCH... files differ. Something might be wrong..."
        exit 2
    else
        echo "OK."
    fi

done

## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm -rf $SFTdir_40h $SFTdir_5d
    rm -f $gridFile Fstats Fstats.log
    for n in 0 1 2 3 6 8 9; do
        rm -f testCFSv2_grid${n}.dat
    done
fi
