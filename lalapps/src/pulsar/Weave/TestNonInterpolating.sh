# Perform a non-interpolating search, and compare F-statistics to reference results

export LAL_FSTAT_FFT_PLAN_MODE=ESTIMATE

echo "=== Create search setup with 3 segments spanning ~3.6 days ==="
set -x
${builddir}/lalapps_WeaveSetup --ref-time=1122334444 --first-segment=1122332211/88000 --segment-count=3 --segment-gap=25000 --detectors=H1,L1 --output-file=WeaveSetup.fits
set +x
echo

echo "=== Restrict timestamps to segment list in WeaveSetup.fits ==="
set -x
${fitsdir}/lalapps_fits_table_list 'WeaveSetup.fits[segments][col c1=start_s; col2=end_s]' \
    | awk 'BEGIN { print "/^#/ { print }" } /^#/ { next } { printf "%i <= $1 && $1 <= %i { print }\n", $1, $2 + 1 }' > timestamp-filter.awk
awk -f timestamp-filter.awk ${srcdir}/timestamps-1.txt > timestamps-1.txt
awk -f timestamp-filter.awk ${srcdir}/timestamps-2.txt > timestamps-2.txt
set +x
echo

echo "=== Extract reference time from WeaveSetup.fits ==="
set -x
${fitsdir}/lalapps_fits_header_getval "WeaveSetup.fits[0]" 'DATE-OBS GPS' > tmp
ref_time=`cat tmp | xargs printf "%.9f"`
set +x
echo

echo "=== Generate SFTs with injected signal ==="
set -x
inject_params="Alpha=2.9; Delta=0.71; Freq=45.5; f1dot=-3e-10"
${injdir}/lalapps_Makefakedata_v5 --randSeed=2345 --fmin=45.0 --Band=1.0 \
    --injectionSources="{refTime=${ref_time}; h0=0.5; cosi=0.1; psi=4.4; phi0=2.1; ${inject_params}}" \
    --Tsft=1800 --outSingleSFT --outSFTdir=. --IFOs=H1,L1 --sqrtSX=1,1 \
    --timestampsFiles=timestamps-2.txt,timestamps-1.txt
set +x
echo

echo "=== Perform non-interpolating search ==="
set -x
${builddir}/lalapps_Weave --output-file=WeaveOut.fits \
    --toplists=mean2F --toplist-limit=2321 --per-detector --per-segment --misc-info \
    --setup-file=WeaveSetup.fits --sft-files='*.sft' \
    --alpha=2.3/0.9 --delta=-1.2/2.3 --freq=45.5~1e-4 --f1dot=-1e-9,0 --semi-max-mismatch=7 --interpolation=no
set +x
echo

echo "=== Check for non-singular semicoherent dimensions ==="
set -x
semi_ntmpl_prev=1
for dim in SSKYA SSKYB NU1DOT NU0DOT; do
    ${fitsdir}/lalapps_fits_header_getval "WeaveOut.fits[0]" "NSEMITMPL ${dim}" > tmp
    semi_ntmpl=`cat tmp | xargs printf "%d"`
    expr ${semi_ntmpl} '>' ${semi_ntmpl_prev}
    semi_ntmpl_prev=${semi_ntmpl}
done
set +x
echo

echo "=== Check computed number of coherent results ==="
set -x
${fitsdir}/lalapps_fits_header_getval "WeaveOut.fits[0]" 'NSEGMENT' > tmp
nsegments=`cat tmp | xargs printf "%d"`
${fitsdir}/lalapps_fits_header_getval "WeaveOut.fits[0]" 'NCOHRES' > tmp
coh_nres=`cat tmp | xargs printf "%d"`
${fitsdir}/lalapps_fits_header_getval "WeaveOut.fits[0]" 'NSEMIRES' > tmp
semi_nres=`cat tmp | xargs printf "%d"`
expr ${coh_nres} '=' ${semi_nres} '*' ${nsegments}
set +x
echo

for seg in 1 2 3; do

    echo "=== Segment #${seg}: Check that no results were recomputed ==="
    set -x
    ${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[per_seg_info][col coh_nrecomp][#row == ${seg}]" > tmp
    coh_nrecomp=`cat tmp | sed "/^#/d" | xargs printf "%d"`
    expr ${coh_nrecomp} '=' 0
    set +x
    echo

done

### Make updating reference results a little easier ###
mkdir TestNonInterpolating.testdir
cp RefSeg1Exact.txt RefSeg2Exact.txt RefSeg3Exact.txt TestNonInterpolating.testdir/
cp WeaveOut.fits TestNonInterpolating.testdir/RefWeaveOut.fits
tar zcvf TestNonInterpolating.tar.gz TestNonInterpolating.testdir/
rm -rf TestNonInterpolating.testdir/

echo "=== Compare semicoherent F-statistics from lalapps_Weave to reference results ==="
set -x
env LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},info" ${builddir}/lalapps_WeaveCompare --setup-file=WeaveSetup.fits --result-file-1=WeaveOut.fits --result-file-2=RefWeaveOut.fits
set +x
echo

echo "=== Compare F-statistic at exact injected signal parameters with loudest F-statistic found by lalapps_Weave ==="
set -x
coh2F_exact=`paste RefSeg1Exact.txt RefSeg2Exact.txt RefSeg3Exact.txt | sed -n '/^[^%]/p' | sed -n '1p' | awk '{print ($7 + $16 + $25) / 3}'`
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean2F_toplist][col c1=mean2F][#row == 1]" > tmp
coh2F_loud=`cat tmp | sed "/^#/d" | xargs printf "%g"`
# Value of 'mean_mu' was calculated by:
#   octapps_run WeaveFstatMismatch --setup-file=TestNonInterpolating.testdir/WeaveSetup.fits --spindowns=1 --semi-max-mismatch=7 --coh-max-mismatch=0 --printarg=meanOfHist
mean_mu=0.53991
awk "BEGIN { print mu = ( ${coh2F_exact} - ${coh2F_loud} ) / ${coh2F_exact}; exit ( mu < ${mean_mu} ? 0 : 1 ) }"
set +x
echo
