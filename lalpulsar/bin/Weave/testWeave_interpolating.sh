if test "${CFITSIO_ENABLED}" = false; then
    echo "Skipping test: requires CFITSIO"
    exit 77
fi
if test "${HOST_OS_DARWIN}" = true; then
    echo "Skipping test; results computed on MacOSX are too different to reference results (computed on Ubuntu Linux)"
    exit 77
fi
case "${CC}" in
    icc*)
        echo "Skipping test; results computed with '${CC}' compiler are too different to reference results (computed with GCC compiler)"
        exit 77
        ;;
    *)
esac

# Perform an interpolating search, and compare F-statistics to reference results

export LAL_FSTAT_FFT_PLAN_MODE=ESTIMATE

echo "=== Create search setup with 3 segments spanning ~260 days ==="
set -x
weavesetup=RefWeaveSetup1.fits
if [ ! -f "${weavesetup}" ]; then
    lalpulsar_WeaveSetup --ephem-earth=earth00-40-DE405.dat.gz --ephem-sun=sun00-40-DE405.dat.gz --first-segment=1122332211/90000 --segment-count=3 --segment-gap=11130000 --detectors=H1,L1 --output-file=${weavesetup}
fi
lalpulsar_fits_overview ${weavesetup}
set +x
echo

echo "=== Restrict timestamps to segment list in ${weavesetup} ==="
set -x
lalpulsar_fits_table_list "${weavesetup}"'[segments][col c1=start_s; col2=end_s]' \
    | awk 'BEGIN { print "/^#/ { print }" } /^#/ { next } { printf "%i <= $1 && $1 <= %i { print }\n", $1, $2 + 1 }' > timestamp-filter.awk
awk -f timestamp-filter.awk all-timestamps-1.txt > timestamps-1.txt
awk -f timestamp-filter.awk all-timestamps-2.txt > timestamps-2.txt
set +x
echo

echo "=== Extract reference time from ${weavesetup} ==="
set -x
ref_time=`lalpulsar_fits_header_getval "${weavesetup}[0]" 'DATE-OBS GPS' | tr '\n\r' '  ' | awk 'NF == 1 {printf "%.9f", $1}'`
set +x
echo

echo "=== Generate SFTs with injected signal ==="
set -x
inject_params="Alpha=2.324; Delta=-1.204; Freq=50.5; f1dot=-1.5e-10"
lalpulsar_Makefakedata_v5 --randSeed=3456 --fmin=50.3 --Band=0.4 --Tsft=1800 \
    --injectionSources="{refTime=${ref_time}; h0=0.5; cosi=0.7; psi=2.1; phi0=3.1; ${inject_params}}" \
    --outSingleSFT --outSFTdir=. --IFOs=H1,L1 --sqrtSX=1,1 \
    --timestampsFiles=timestamps-1.txt,timestamps-2.txt
set +x
echo

echo "=== Perform interpolating search ==="
set -x
lalpulsar_Weave --output-file=WeaveOut.fits \
    --toplists=mean2F --toplist-limit=2321 \
    --extra-statistics="coh2F,coh2F_det,mean2F_det" \
    --segment-info --time-search --setup-file=${weavesetup} --sft-files='*.sft' \
    --Fstat-method=ResampBest --Fstat-run-med-window=50 \
    --alpha=2.3/0.05 --delta=-1.2/0.1 --freq=50.5~0.005 --f1dot=-3e-10,0 \
    --semi-max-mismatch=6.5 --coh-max-mismatch=0.4
lalpulsar_fits_overview WeaveOut.fits
set +x
echo

echo "=== Check for non-singular semicoherent dimensions ==="
set -x
semi_ntmpl_prev=1
for dim in SSKYA SSKYB NU1DOT NU0DOT; do
    semi_ntmpl=`lalpulsar_fits_header_getval "WeaveOut.fits[0]" "NSEMITMPL ${dim}" | tr '\n\r' '  ' | awk 'NF == 1 {printf "%d", $1}'`
    expr ${semi_ntmpl} '>' ${semi_ntmpl_prev}
    semi_ntmpl_prev=${semi_ntmpl}
done
set +x
echo

weave_recomp_threshold=0.0
echo "=== Check that number of recomputed results is below tolerance ==="
set -x
coh_nres=`lalpulsar_fits_header_getval "WeaveOut.fits[0]" 'NCOHRES' | tr '\n\r' '  ' | awk 'NF == 1 {printf "%d", $1}'`
coh_ntmpl=`lalpulsar_fits_header_getval "WeaveOut.fits[0]" 'NCOHTPL' | tr '\n\r' '  ' | awk 'NF == 1 {printf "%d", $1}'`
awk "BEGIN { print recomp = ( ${coh_nres} - ${coh_ntmpl} ) / ${coh_ntmpl}; exit ( recomp <= ${weave_recomp_threshold} ? 0 : 1 ) }"
expr ${coh_nres} '=' ${coh_ntmpl}
set +x
echo

### Make updating reference results a little easier ###
mkdir newtarball/
cd newtarball/
cp ../${weavesetup} .
cp ../all-timestamps-1.txt ../all-timestamps-2.txt .
cp ../RefSeg1Exact.txt ../RefSeg2Exact.txt ../RefSeg3Exact.txt .
cp ../WeaveOut.fits RefWeaveOut.fits
tar zcf ../new_testWeave_interpolating.tar.gz *
cd ..
rm -rf newtarball/

echo "=== Compare semicoherent F-statistics from lalpulsar_Weave to reference results ==="
set -x
env LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},info" lalpulsar_WeaveCompare --setup-file=${weavesetup} --result-file-1=WeaveOut.fits --result-file-2=RefWeaveOut.fits \
    --round-param-to-dp=9 --round-param-to-sf=3 --toplist-limit=25 --result-tol-L2=8e-2 --result-tol-angle=8e-2
set +x
echo

echo "=== Compare F-statistic at exact injected signal parameters with loudest F-statistic found by lalpulsar_Weave ==="
set -x
coh2F_exact=`paste RefSeg1Exact.txt RefSeg2Exact.txt RefSeg3Exact.txt | sed -n '/^[^%]/p' | sed -n '1p' | awk '{print ($7 + $16 + $25) / 3}'`
lalpulsar_fits_table_list "WeaveOut.fits[mean2F_toplist][col c1=mean2F][#row == 1]" > tmp
coh2F_loud=`cat tmp | sed "/^#/d" | xargs printf "%.16g"`
# Value of 'mean_mu' was calculated by:
#   octapps_run WeaveFstatMismatch --setup-file=${weavesetup} --spindowns=1 --semi-max-mismatch=5.5 --coh-max-mismatch=0.3 --printarg=meanOfHist
mean_mu=0.47315
awk "BEGIN { print mu = ( ${coh2F_exact} - ${coh2F_loud} ) / ${coh2F_exact}; exit ( mu < ${mean_mu} ? 0 : 1 ) }"
set +x
echo

echo "=== Check that lalpulsar_Weave can be run at a single point ==="
set -x
lalpulsar_fits_table_list "WeaveOut.fits[mean2F_toplist][col c1=alpha; c2=delta; c3=freq; c4=f1dot; c5=mean2F][#row == 1]" > tmp
alpha_loud=`cat tmp | sed "/^#/d" | awk '{print $1}' | xargs printf "%.16g"`
delta_loud=`cat tmp | sed "/^#/d" | awk '{print $2}' | xargs printf "%.16g"`
freq_loud=` cat tmp | sed "/^#/d" | awk '{print $3}' | xargs printf "%.16g"`
f1dot_loud=`cat tmp | sed "/^#/d" | awk '{print $4}' | xargs printf "%.16g"`
coh2F_loud=`cat tmp | sed "/^#/d" | awk '{print $5}' | xargs printf "%.16g"`
lalpulsar_Weave --output-file=WeaveOutSingle.fits \
    --toplists=mean2F,log10BSGL,log10BSGLtL,log10BtSGLtL --toplist-limit=2321 \
    --extra-statistics="coh2F,coh2F_det,mean2F_det,ncount,ncount_det" --lrs-Fstar0sc=2000 --lrs-oLGX=4,0.1 \
    --setup-file=${weavesetup} --sft-files='*.sft' \
    --alpha=${alpha_loud}~0 --delta=${delta_loud}~0 --freq=${freq_loud}~0 --f1dot=${f1dot_loud}~0 \
    --semi-max-mismatch=6.5 --coh-max-mismatch=0.4
lalpulsar_fits_overview WeaveOutSingle.fits
lalpulsar_fits_table_list "WeaveOutSingle.fits[mean2F_toplist][col c1=mean2F][#row == 1]" > tmp
coh2F_loud_single=`cat tmp | sed "/^#/d" | xargs printf "%.16g"`
mean_mu=0.15
awk "BEGIN { print mu = ( ${coh2F_loud} - ${coh2F_loud_single} ) / ${coh2F_loud}; exit ( mu < ${mean_mu} ? 0 : 1 ) }"
set +x
echo
