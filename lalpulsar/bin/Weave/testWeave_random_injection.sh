if test "${CFITSIO_ENABLED}" = false; then
    echo "Skipping test: requires CFITSIO"
    exit 77
fi

# Perform an interpolating search with a random injection, and compare to injection generated with lalpulsar_Makefakedata_v5

export LAL_FSTAT_FFT_PLAN_MODE=ESTIMATE

echo "=== Create search setup with 3 segments spanning ~3.6 days ==="
set -x
lalpulsar_WeaveSetup --ephem-earth=earth00-19-DE405.dat.gz --ephem-sun=sun00-19-DE405.dat.gz --ref-time=1122334444 --first-segment=1122332211/88000 --segment-count=3 --segment-gap=25000 --detectors=H1,L1 --output-file=WeaveSetup.fits
lalpulsar_fits_overview WeaveSetup.fits
set +x
echo

echo "=== Restrict timestamps to segment list in WeaveSetup.fits ==="
set -x
lalpulsar_fits_table_list 'WeaveSetup.fits[segments][col c1=start_s; col2=end_s]' \
    | awk 'BEGIN { print "/^#/ { print }" } /^#/ { next } { printf "%i <= $1 && $1 <= %i { print }\n", $1, $2 + 1 }' > timestamp-filter.awk
awk -f timestamp-filter.awk all-timestamps-1.txt > timestamps-1.txt
awk -f timestamp-filter.awk all-timestamps-2.txt > timestamps-2.txt
set +x
echo

echo "=== Extract reference time from WeaveSetup.fits ==="
set -x
ref_time=`lalpulsar_fits_header_getval "WeaveSetup.fits[0]" 'DATE-OBS GPS' | tr '\n\r' '  ' | awk 'NF == 1 {printf "%.9f", $1}'`
set +x
echo

echo "=== Generate SFTs with no injected signal ==="
set -x
mkdir -p no_inj/
lalpulsar_Makefakedata_v5 --randSeed=3456 --fmin=50.3 --Band=0.4 --Tsft=1800 \
    --outSingleSFT --outSFTdir=no_inj/ --IFOs=H1,L1 --sqrtSX=1,1 \
    --timestampsFiles=timestamps-1.txt,timestamps-2.txt
set +x
echo

echo "=== Perform interpolating search with random injection ==="
set -x
lalpulsar_Weave --rand-seed=9876 --output-file=WeaveOutRandInj.fits \
    --random-injection=0.5 \
    --toplists=mean2F --toplist-limit=10 \
    --extra-statistics="coh2F,coh2F_det,mean2F_det" \
    --segment-info --time-search \
    --setup-file=WeaveSetup.fits --sft-files='no_inj/*.sft' \
    --Fstat-method=ResampBest --Fstat-run-med-window=50 \
    --alpha=2.3/0.05 --delta=-1.2/0.1 --freq=50.5~0.005 --f1dot=-3e-10,0 \
    --semi-max-mismatch=6.5 --coh-max-mismatch=0.4
lalpulsar_fits_overview WeaveOutRandInj.fits
set +x
echo

echo "=== Extract random injection ==="
set -x
echo "[random_injection]" > random_injection.cfg
echo "refTime = 1122334444" >> random_injection.cfg
for par in psi phi0 aPlus aCross Alpha Delta Freq f1dot; do
    lalpulsar_fits_table_list "WeaveOutRandInj.fits[injections][col c1=${par}][#row == 1]" > tmp
    cat tmp | sed "/^#/d" | awk '{print $1}' | xargs printf "${par} = %.16g\n" >> random_injection.cfg
done
cat random_injection.cfg
set +x
echo

echo "=== Add extracted injected signal to SFTs ==="
set -x
mkdir -p inj/
lalpulsar_Makefakedata_v5  --noiseSFTs='no_inj/*.sft' --SFTWindowType=rectangular \
    --injectionSources=./random_injection.cfg \
    --outSingleSFT --outSFTdir=inj/
set +x
echo

echo "=== Perform interpolating search with Makefakedata_v5 injection ==="
set -x
lalpulsar_Weave --output-file=WeaveOutMFDInj.fits \
    --toplists=mean2F --toplist-limit=10 \
    --extra-statistics="coh2F,coh2F_det,mean2F_det" \
    --segment-info --time-search \
    --setup-file=WeaveSetup.fits --sft-files='inj/*.sft' \
    --Fstat-method=ResampBest --Fstat-run-med-window=50 \
    --alpha=2.3/0.05 --delta=-1.2/0.1 --freq=50.5~0.005 --f1dot=-3e-10,0 \
    --semi-max-mismatch=6.5 --coh-max-mismatch=0.4
lalpulsar_fits_overview WeaveOutMFDInj.fits
set +x
echo

echo "=== Compare lalpulsar_Weave results ==="
set -x
LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},info" lalpulsar_WeaveCompare --setup-file=WeaveSetup.fits --result-file-1=WeaveOutRandInj.fits --result-file-2=WeaveOutMFDInj.fits --sort-by-semi-phys
set +x
echo
