#!/bin/sh


## take user-arguments
extra_args="$@"

## allow 'make test' to work from builddir != srcdir
if [ -n "${srcdir}" ]; then
    builddir="./";
else
    srcdir=.
fi

##---------- names of codes and input/output files
psd_code="${builddir}lalapps_ComputePSD"
SFTdir="$srcdir"

# test if LAL_DATA_PATH has been set ... needed to locate ephemeris-files
if [ -z "$LAL_DATA_PATH" ]; then
    if [ -n "$LALPULSAR_PREFIX" ]; then
	export LAL_DATA_PATH=".:${LALPULSAR_PREFIX}/share/lalpulsar";
    else
	echo
	echo "Need environment-variable LALPULSAR_PREFIX, or LAL_DATA_PATH to be set"
	echo "to your ephemeris-directory (e.g. /usr/local/share/lalpulsar)"
	echo "This might indicate an incomplete LAL+LALPULSAR installation"
	echo
	exit 1
    fi
fi

tolerance=1e-5

# ---------- fixed parameter of our test-signal
IFO=H1

fStart=300.0
fBand=0.1


blocksRngMed=101
outPSD=psd1.dat
refPSD=${srcdir}/psd_ref.dat

inputData="${SFTdir}/SFT.0*"


## ----- run computePSD
cmdline="${psd_code} --IFO=$IFO --inputData='$inputData' --outputPSD=$outPSD --blocksRngMed=$blocksRngMed --fStart=$fStart --fBand=$fBand"

echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$psd_code' ..."
    exit 1
fi

## ----- compare result PSD to reference result in source-directory
if [ ! -r "$outPSD" -o ! -r $ "$refPSD" ]; then
    echo "ERROR: missing psd output file '$outPSD' or '$refPSD'"
    exit 1
fi

cmp=`paste $outPSD $refPSD | LC_ALL=C awk 'BEGIN {n=0; maxErr = 0; avgErr = 0; binsOff = 0} {n+=1; dFreq = $3 - $1; if ( dFreq != 0 ) binsOff ++; relErr = 2*($4 - $2)/($4 + $2); if (relErr > maxErr) maxErr = relErr; avgErr += relErr } END { avgErr /= n; printf "binsOff=%d; avgErr=%g; maxErr=%g", binsOff, avgErr, maxErr}'`

eval $cmp

failed=`echo $maxErr $tolerance | awk '{failed = $1 > $2; print failed}'`

if [ "$binsOff" != "0" -o "$failed" != "0" ]; then
    echo
    echo "*****   ComputePSD Test FAILED *****"
    echo "binsOff = $binsOff (tolerance = 0)"
    echo "avgErr = $avgErr"
    echo "maxErr = $maxErr  (tolerance = $tolerance)"
    echo "(maxErr > tolerance) = $failed"
    echo "************************************"
    echo
    exit 1
else
    echo
    echo "========== OK. ComputePSD test PASSSED. =========="
    echo
    exit 0
fi
