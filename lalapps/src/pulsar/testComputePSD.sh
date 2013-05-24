#!/bin/sh


## take user-arguments
extra_args="$@"

## allow 'make test' to work from builddir != srcdir
if [ -z "${srcdir}" ]; then
    srcdir=`dirname $0`
fi

builddir="./";
injectdir="./Injections/"

## check for LALPulsar data directory
if [ -z "${LALPULSAR_DATADIR}" ]; then
    echo "Need environment-variable LALPULSAR_DATADIR to be set to" 
    echo "your ephemeris-directory (e.g. /usr/local/share/lalpulsar)" 
    echo "This might indicate an incomplete LAL+LALPULSAR installation" 
    exit 1
fi

##---------- names of codes and input/output files
psd_code="${builddir}lalapps_ComputePSD"
mfd_code="${injectdir}lalapps_Makefakedata_v4 -E ${LALPULSAR_DATADIR}" 
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
outPSDstripped=psd1-stripped.dat
refPSD=${srcdir}/psd_ref.dat

inputData="${SFTdir}/SFT.0*"

echo "ComputePSD Test1: PSD comparison to reference result..."

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
cat ${outPSD} | sed -e"/^%%.*/d" > ${outPSDstripped}

cmp=`paste ${outPSDstripped} $refPSD | LC_ALL=C awk 'BEGIN {n=0; maxErr = 0; avgErr = 0; binsOff = 0} {n+=1; dFreq = $3 - $1; if ( dFreq != 0 ) binsOff ++; relErr = 2*($4 - $2)/($4 + $2); if (relErr > maxErr) maxErr = relErr; avgErr += relErr } END { avgErr /= n; printf "binsOff=%d; avgErr=%g; maxErr=%g", binsOff, avgErr, maxErr}'`

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
    echo "***** ComputePSD Test1 passed *****"
    echo
fi

## ----- Test 2 - correct cropping of normalized SFT

echo "ComputePSD Test2: --outputNormSFT frequency range..."

outSFT="./testpsd_sft_H1"
linefreq="50.05"

## ----- run MFDv4
cmdline="${mfd_code} --IFO=$IFO --outSingleSFT=1 --outSFTbname=$outSFT --startTime=828002611 --duration=1800 --fmin=50 --Band=0.1 --noiseSqrtSh=3.25e-22 --lineFeature=1 --h0=5e-23 --cosi=0 --Freq=$linefreq --randSeed=1 "

echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

outPSD_full="./psd_fullsft.dat"
outPSD_band="./psd_freqband.dat"

## ----- run computePSD once for full SFT
cmdline="${psd_code} --inputData=$outSFT --outputPSD=$outPSD_full --blocksRngMed=$blocksRngMed --outputNormSFT=1"

echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$psd_code' ..."
    exit 1
fi

## ----- run computePSD for a second time, with restricted band
cmdline="${psd_code} --inputData=$outSFT --outputPSD=$outPSD_band --blocksRngMed=$blocksRngMed --outputNormSFT=1  --Freq=50.03 --FreqBand=0.04"

echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$psd_code' ..."
    exit 1
fi

topline_psd_full=$(sort -nr -k2,2 $outPSD_full | head -1)
toppsd_full=$(echo $topline_psd_full | awk '{print $2}')
toppsdfreq_full=$(echo $topline_psd_full | awk '{print $1}')
topline_power_full=$(sort -nr -k3,3 $outPSD_full | head -1)
toppower_full=$(echo $topline_power_full | awk '{print $3}')
toppowerfreq_full=$(echo $topline_power_full | awk '{print $1}')

topline_psd_band=$(sort -nr -k2,2 $outPSD_band | head -1)
toppsd_band=$(echo $topline_psd_band | awk '{print $2}')
toppsdfreq_band=$(echo $topline_psd_band | awk '{print $1}')
topline_power_band=$(sort -nr -k3,3 $outPSD_band | head -1)
toppower_band=$(echo $topline_power_band | awk '{print $3}')
toppowerfreq_band=$(echo $topline_power_band | awk '{print $1}')

echo "Loudest bins:"      
echo "==>  full SFT: PSD=$toppsd_full at $toppsdfreq_full Hz, normPower=$toppower_full at $toppowerfreq_full Hz"
echo "==>  freqband: PSD=$toppsd_band at $toppsdfreq_band Hz, normPower=$toppower_band at $toppowerfreq_band Hz"

retstatus=0
tolerance_freq=0.0
tolerance_psd=1e-46
tolerance_power=0.001
tolerance_rel=0.001
awk_absdiff='{printf "%.6f", sqrt(($1-$2)*($1-$2)) }'
awk_reldev='{printf "%.2e", sqrt(($1-$2)*($1-$2))/(0.5*($1+$2)) }'
awk_isgtr='{if($1>$2) {print "1"}}'
freqdiff_full=$(echo $toppowerfreq_full $linefreq | awk "$awk_absdiff")
freqdiff_band=$(echo $toppowerfreq_band $linefreq | awk "$awk_absdiff")
freqdiff_runs=$(echo $toppowerfreq_full $toppowerfreq_band | awk "$awk_absdiff")
toppsd_reldev=$(echo $toppsd_full $toppsd_band | awk "$awk_reldev")
toppower_reldev=$(echo $toppower_full $toppower_band | awk "$awk_reldev")
echo "Injected line recovered in full run with offset $freqdiff_full Hz and in band run with $freqdiff_band Hz - difference between runs: $freqdiff_runs Hz."
echo "Relative difference in highest PSD value is $toppsd_reldev and in highest normPower value $toppower_reldev."

fail1=$(echo $freqdiff_full $tolerance_freq | awk "$awk_isgtr")
fail2=$(echo $freqdiff_band $tolerance_freq | awk "$awk_isgtr")
fail3=$(echo $freqdiff_runs $tolerance_freq | awk "$awk_isgtr")
fail4=$(echo $toppsd_reldev $tolerance_rel  | awk "$awk_isgtr")
fail5=$(echo $toppower_reldev $tolerance_rel  | awk "$awk_isgtr")
if [ "$fail1" -o "$fail2" -o "$fail3" -o "$fail4" -o "$fail5" ]; then
    echo " ==> FAILED"
    retstatus=1
else
    echo " ==> OK"
    echo
    echo "========== OK. All ComputePSD tests PASSED. =========="
    echo
fi

## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm $outPSD $outPSDstripped
    rm $outSFT
    rm $outPSD_band
    rm $outPSD_full
    echo "Cleaned up."
fi

exit $retstatus
