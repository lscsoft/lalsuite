#!/bin/sh

## allow 'make test' to work from builddir != srcdir
if [ -n "${srcdir}" ]; then
    builddir="./";
else
    srcdir=.
fi

sftdir="${srcdir}/.."

sftbase="SFT.0000"
IFO="LHO"
FCOMPARE="${builddir}lalapps_compareFstats"
CFS_DEFAULT="${builddir}lalapps_ComputeFStatistic"

outfile1="Fstatv1_1.dat";
outfile2="Fstatv1_2.dat";
outfile3="Fstatv1_3.dat";

CFSparams1="--IFO=$IFO --DataDir=$sftdir --BaseName=$sftbase --Freq=300.1 --Fthreshold=0\
--FreqBand=0.2 --Alpha=2.2 --AlphaBand=0.012 --Delta=0.8 --DeltaBand=0.018 --gridType=0 --outputFstat=$outfile1"

CFSparams2="--IFO=$IFO --DataDir=$sftdir --BaseName=$sftbase --Freq=300.1 --Fthreshold=0\
--FreqBand=0.2 --Alpha=2.2 --AlphaBand=0.003 --Delta=0.8 --DeltaBand=0.003 --gridType=1 --outputFstat=$outfile2"

#give help string if requested
if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
    echo
    echo "Usage: $0 [yourCFScode]"
    echo
    echo "The default-code used is '$CFS_DEFAULT'"
    echo
    exit 1
fi

if [ x$1 = x ]; then
    prog=$CFS_DEFAULT;
    extra_args=
else
    prog=$1;
    shift
    extra_args="$@"
fi


# test if LAL_DATA_PATH has been set ... needed to locate ephemeris-files
if [ -z "$LAL_DATA_PATH" ]; then
    echo
    echo "Need environment-variable LAL_DATA_PATH to point to your ephemeris-directory (e.g. /usr/local/share/lal)"
    if [ -n "$LAL_PREFIX" ]; then
	echo "You have LAL_PREFIX set, I suggest setting 'LAL_DATA_PATH=\$LAL_PREFIX/share/lal'"
    fi
    echo
    exit 1
fi

## Tests start here
## --------------------
echo
echo "Running ComputeFStatistic-code '$prog' on test-data '$sftdir/$sftbase*'"

## Test1: using a uniform sky-grid
##----------------------------------------
echo
echo "----------------------------------------------------------------------"
echo "Test 1) uniform sky-grid:"
echo "----------------------------------------------------------------------"

cmdline="$prog $CFSparams1 $extra_args";
echo $cmdline
if ! $cmdline ; then
    echo "Something failed ... giving up.";
   exit 2;
fi

echo
echo "Comparing output-file 'Fstats' with reference-version 'Fstats.ref1' ... "

cmdline="$FCOMPARE --clusterFiles=false -1 ./${outfile1} -2 ${srcdir}/Fstats.ref1 --Ftolerance=0.01";
echo $cmdline
if $cmdline &> test1.dat; then
    echo "OK."
else
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
fi

## Test2: using an isotropic Grid
##-------------------------------
echo
echo "----------------------------------------------------------------------"
echo "Test 2) isotropic sky-grid:"
echo "----------------------------------------------------------------------"

cmdline="$prog $CFSparams2 $extra_args"
echo $cmdline
if ! $cmdline; then
    echo "Something failed ... giving up.";
    exit 2;
fi

echo
echo "Comparing output-file 'Fstats' with reference-version 'Fstats.ref2' ... "
cmdline="$FCOMPARE --clusterFiles=false -1 ./${outfile2} -2 ${srcdir}/Fstats.ref2 --Ftolerance=0.01"
echo $cmdline
if $cmdline &> test2.dat; then
    echo "OK."
else
    echo "OUCH... files differ. Something might be wrong..."
    exit 2
fi
