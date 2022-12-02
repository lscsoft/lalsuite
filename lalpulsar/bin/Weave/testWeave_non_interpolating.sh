if test "${CFITSIO_ENABLED}" = false; then
    echo "Skipping test: requires CFITSIO"
    exit 77
fi

# Perform a non-interpolating search, and compare F-statistics to reference results

export LAL_FSTAT_FFT_PLAN_MODE=ESTIMATE

echo "=== Create search setup with 3 segments spanning ~3.6 days ==="
set -x
weavesetup=RefWeaveSetup1.fits
if [ ! -f "${weavesetup}" ]; then
    lalpulsar_WeaveSetup --ephem-earth=earth00-40-DE405.dat.gz --ephem-sun=sun00-40-DE405.dat.gz --ref-time=1122334444 --first-segment=1122332211/88000 --segment-count=3 --segment-gap=25000 --detectors=H1,L1 --output-file=${weavesetup}
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
inject_params="Alpha=2.9; Delta=0.71; Freq=45.5; f1dot=-3e-10"
lalpulsar_Makefakedata_v5 --randSeed=2345 --fmin=45.0 --Band=1.0 \
    --injectionSources="{refTime=${ref_time}; h0=0.5; cosi=0.1; psi=4.4; phi0=2.1; ${inject_params}}" \
    --Tsft=1800 --outSingleSFT --outSFTdir=. --IFOs=H1,L1 --sqrtSX=1,1 \
    --timestampsFiles=timestamps-2.txt,timestamps-1.txt
set +x
echo

echo "=== Perform non-interpolating search ==="
set -x
lalpulsar_Weave --output-file=WeaveOut.fits \
    --toplists=mean2F --toplist-limit=2321 --extra-statistics="mean2F_det,coh2F,coh2F_det" \
    --segment-info --setup-file=${weavesetup} --sft-files='*.sft' \
    --Fstat-method=ResampGeneric --Fstat-run-med-window=50 \
    --alpha=2.3/0.9 --delta=-1.2/2.3 --freq=45.5~0.005 --f1dot=-1e-9,0 \
    --semi-max-mismatch=7 --interpolation=no
lalpulsar_fits_overview WeaveOut.fits
set +x
echo

echo "=== Perform non-interpolating search with recalc ==="
set -x
lalpulsar_Weave --output-file=WeaveOutRecalc.fits \
    --toplists=mean2F --toplist-limit=2321 --extra-statistics=all --recalc-statistics=all --Fstat-method=DemodBest \
    --setup-file=${weavesetup} --sft-files='*.sft' \
    --Fstat-run-med-window=50 \
    --alpha=2.3/0.9 --delta=-1.2/2.3 --freq=45.5~0.005 --f1dot=-1e-9,0 \
    --semi-max-mismatch=7 --interpolation=no
lalpulsar_fits_overview WeaveOutRecalc.fits
set +x
echo

echo "=== Check recalc'ed results against original ==="
set -x
lalpulsar_fits_table_list -n "WeaveOutRecalc.fits[mean2F_toplist][col coh2F;     coh2F_H1;     coh2F_L1; ]"   > stage0.dat
lalpulsar_fits_table_list -n "WeaveOutRecalc.fits[mean2F_toplist][col coh2F_rec; coh2F_H1_rec; coh2F_L1_rec]" > stage1.dat
diff stage0.dat stage1.dat
set +x

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

echo "=== Check that no results were recomputed ==="
set -x
coh_nres=`lalpulsar_fits_header_getval "WeaveOut.fits[0]" 'NCOHRES' | tr '\n\r' '  ' | awk 'NF == 1 {printf "%d", $1}'`
coh_ntmpl=`lalpulsar_fits_header_getval "WeaveOut.fits[0]" 'NCOHTPL' | tr '\n\r' '  ' | awk 'NF == 1 {printf "%d", $1}'`
expr ${coh_nres} '=' ${coh_ntmpl}
set +x
echo

echo "=== Check computed number of coherent and semicoherent templates ==="
set -x
nsegments=`lalpulsar_fits_header_getval "WeaveOut.fits[0]" 'NSEGMENT' | tr '\n\r' '  ' | awk 'NF == 1 {printf "%d", $1}'`
coh_ntmpl=`lalpulsar_fits_header_getval "WeaveOut.fits[0]" 'NCOHTPL' | tr '\n\r' '  ' | awk 'NF == 1 {printf "%d", $1}'`
semi_ntmpl=`lalpulsar_fits_header_getval "WeaveOut.fits[0]" 'NSEMITPL' | tr '\n\r' '  ' | awk 'NF == 1 {printf "%d", $1}'`
expr ${coh_ntmpl} '=' ${semi_ntmpl} '*' ${nsegments}
set +x
echo

### Make updating reference results a little easier ###
mkdir newtarball/
cd newtarball/
cp ../${weavesetup} .
cp ../all-timestamps-1.txt ../all-timestamps-2.txt .
cp ../RefSeg1Exact.txt ../RefSeg2Exact.txt ../RefSeg3Exact.txt .
cp ../WeaveOut.fits RefWeaveOut.fits
tar zcf ../new_testWeave_non_interpolating.tar.gz *
cd ..
rm -rf newtarball/

echo "=== Compare semicoherent F-statistics from lalpulsar_Weave to reference results ==="
set -x
env LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},info" lalpulsar_WeaveCompare --setup-file=${weavesetup} --result-file-1=WeaveOut.fits --result-file-2=RefWeaveOut.fits \
    --round-param-to-dp=10 --round-param-to-sf=4 --unmatched-item-tol=0.5
set +x
echo

echo "=== Compare F-statistic at exact injected signal parameters with loudest F-statistic found by lalpulsar_Weave ==="
set -x
coh2F_exact=`paste RefSeg1Exact.txt RefSeg2Exact.txt RefSeg3Exact.txt | sed -n '/^[^%]/p' | sed -n '1p' | awk '{print ($7 + $16 + $25) / 3}'`
lalpulsar_fits_table_list "WeaveOut.fits[mean2F_toplist][col c1=mean2F][#row == 1]" > tmp
coh2F_loud=`cat tmp | sed "/^#/d" | xargs printf "%.16g"`
# Value of 'mean_mu' was calculated by:
#   octapps_run WeaveFstatMismatch --setup-file=${weavesetup} --spindowns=1 --semi-max-mismatch=7 --coh-max-mismatch=0 --printarg=meanOfHist
mean_mu=0.53991
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
    --toplists=mean2F --toplist-limit=2321 --extra-statistics="mean2F_det,coh2F,coh2F_det" \
    --setup-file=${weavesetup} --sft-files='*.sft' \
    --alpha=${alpha_loud}~0 --delta=${delta_loud}~0 --freq=${freq_loud}~0 --f1dot=${f1dot_loud}~0 \
    --semi-max-mismatch=7 --interpolation=no
lalpulsar_fits_overview WeaveOutSingle.fits
lalpulsar_fits_table_list "WeaveOutSingle.fits[mean2F_toplist][col c1=mean2F][#row == 1]" > tmp
coh2F_loud_single=`cat tmp | sed "/^#/d" | xargs printf "%.16g"`
mean_mu=0.05
awk "BEGIN { print mu = ( ${coh2F_loud} - ${coh2F_loud_single} ) / ${coh2F_loud}; exit ( mu < ${mean_mu} ? 0 : 1 ) }"
set +x
echo
