#!/bin/sh
sftdir=".."
sftbase="SFT.0000"
IFO="LHO"

CFSparams="--IFO=$IFO --DataDir=$sftdir --BaseName=$sftbase --Freq=300.1 \
--FreqBand=0.2 --Alpha=2.2 --AlphaBand=0.003 --Delta=0.8 --DeltaBand=0.003"

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

echo "Running ComputeFStatistic-code '$prog' on test-data '$sftdir/$sftbase*'"
if "$prog" $CFSparams -v0; then
    echo "done";
else
    echo "failed... exiting.";
    exit
fi

echo "Comparing output-file 'Fstats' with reference-version 'Fstats.ref' ..."

if ./compareFstats -1 Fstats -2 Fstats.ref ; then
    echo "OK. No differences found!"
else
    echo "OUCH... files differ. Something might be wrong..."
fi

    

