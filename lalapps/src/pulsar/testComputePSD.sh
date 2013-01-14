#!/bin/sh


## take user-arguments
extra_args="$@"

## allow 'make test' to work from builddir != srcdir
if [ -z "${srcdir}" ]; then
    srcdir=`dirname $0`
fi

builddir="./";

##---------- names of codes and input/output files
psd_code="${builddir}lalapps_ComputePSD"
SFTdir="$srcdir"

tolerance=1e-5

# ---------- fixed parameter of our test-signal
IFO=H1

## we specify the frequency range as an "open" interval (300, 300.1) to avoid
## roundoff ambiguities to discrete bins on different platforms, using eps=1e-6
fStart=300.000001 ## 300 + eps
fBand=0.099997	  ## 0.1 - 3eps


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
if [ ! -r "$outPSD" -o ! -r "$refPSD" ]; then
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
