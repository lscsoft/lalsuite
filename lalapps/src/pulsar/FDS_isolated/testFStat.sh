#!/bin/sh
sftdir=".."
sftbase="SFT.0000"
IFO="LHO"

CFSparams1="--IFO=$IFO --DataDir=$sftdir --BaseName=$sftbase --Freq=300.1 \
--FreqBand=0.2 --Alpha=2.2 --AlphaBand=0.003 --Delta=0.8 --DeltaBand=0.003 --gridType=0"

CFSparams2="--IFO=$IFO --DataDir=$sftdir --BaseName=$sftbase --Freq=300.1 \
--FreqBand=0.2 --Alpha=2.2 --AlphaBand=0.003 --Delta=0.8 --DeltaBand=0.003 --gridType=1"

CFSparams3="--IFO=$IFO --DataDir=$sftdir --BaseName=$sftbase --Freq=300.1 \
--FreqBand=0.2 --Alpha=2.2 --AlphaBand=1.0 --Delta=0.8 --DeltaBand=1.0 \
--gridType=2 --metricType=1 --metricMismatch=0.02"

# test if LAL_DATA_PATH has been set ... needed to locate ephemeris-files
if [ x$LAL_DATA_PATH = x ]; then
    echo
    echo "Need environment-variable LAL_DATA_PATH to point to your ephemeris-directory (e.g. /usr/local/share/lal)"
    echo
    exit
fi


if [ x$1 = x ]; then
    prog="./lalapps_ComputeFStatistic";
else
    prog="$1";
fi

## Tests start here 
## --------------------
echo "Running ComputeFStatistic-code '$prog' on test-data '$sftdir/$sftbase*'"

## Test1: using a uniform sky-grid
##----------------------------------------
echo "Test 1) uniform sky-grid:"
echo "cmd-line args: $CFSparams1"
if ! "$prog" $CFSparams1 -v0; then
    echo "failed... exiting.";
    exit
fi

echo -n "Comparing output-file 'Fstats' with reference-version 'Fstats.ref1' ... "

if ./compareFstats -1 ./Fstats -2 ./Fstats.ref1 ; then
    echo "OK."
else
    echo "OUCH... files differ. Something might be wrong..."
fi

## Test2: using an isotropic Grid
##-------------------------------
echo "Test 2) isotropic sky-grid:"
echo "cmd-line args: $CFSparams2"
if ! "$prog" $CFSparams2 -v0; then
    echo "failed... exiting.";
    exit
fi

echo -n "Comparing output-file 'Fstats' with reference-version 'Fstats.ref2' ... "

if ./compareFstats -1 ./Fstats -2 ./Fstats.ref2 ; then
    echo "OK."
else
    echo "OUCH... files differ. Something might be wrong..."
fi


## Test3: using a the analytic Ptole-metric
##----------------------------------------
echo "Test 3) analytic Ptole-metric:"
echo "cmd-line args: $CFSparams3"
if ! "$prog" $CFSparams3 -v0; then
    echo "failed... exiting.";
    exit
fi

echo -n "Comparing output-file 'Fstats' with reference-version 'Fstats.ref3' ... "

if ./compareFstats -1 ./Fstats -2 ./Fstats.ref3 ; then
    echo "OK."
else
    echo "OUCH... files differ. Something might be wrong..."
fi
