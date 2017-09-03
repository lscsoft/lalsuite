# Perform an interpolating search, and compare F-statistics to reference results

export LAL_FSTAT_FFT_PLAN_MODE=ESTIMATE

echo "=== Create search setup with 3 segments spanning ~260 days ==="
set -x
${builddir}/lalapps_WeaveSetup --first-segment=1122332211/90000 --segment-count=3 --segment-gap=11130000 --detectors=H1,L1 --output-file=WeaveSetup.fits
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
inject_params="Alpha=2.324; Delta=-1.204; Freq=50.5; f1dot=-1.5e-10"
${injdir}/lalapps_Makefakedata_v5 --randSeed=3456 --fmin=50.3 --Band=0.4 --Tsft=1800 \
    --injectionSources="{refTime=${ref_time}; h0=0.5; cosi=0.7; psi=2.1; phi0=3.1; ${inject_params}}" \
    --outSingleSFT --outSFTdir=. --IFOs=H1,L1 --sqrtSX=1,1 \
    --timestampsFiles=timestamps-1.txt,timestamps-2.txt
set +x
echo

echo "=== Perform interpolating search ==="
set -x
${builddir}/lalapps_Weave --output-file=WeaveOut.fits \
    --toplists=mean2F,B_S/GL,B_S/GLtL,B_tS/GLtL --toplist-limit=2321 --extra-statistics="coh2F,coh2F_det,mean2F_det,ncount,ncount_det" --lrs-Fstar0sc=2000 --lrs-oLGX=4,0.1 --misc-info \
    --setup-file=WeaveSetup.fits --sft-files='*.sft' \
    --alpha=2.3/0.05 --delta=-1.2/0.1 --freq=50.5~5e-6 --f1dot=-3e-10,0 --semi-max-mismatch=5.5 --coh-max-mismatch=0.3
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

echo "=== Check that no results were recomputed ==="
set -x
${fitsdir}/lalapps_fits_header_getval "WeaveOut.fits[0]" 'NCOHRES' > tmp
coh_nres=`cat tmp | xargs printf "%d"`
${fitsdir}/lalapps_fits_header_getval "WeaveOut.fits[0]" 'NCOHTPL' > tmp
coh_ntmpl=`cat tmp | xargs printf "%d"`
expr ${coh_nres} '=' ${coh_ntmpl}
set +x
echo

### Make updating reference results a little easier ###
mkdir TestInterpolating.testdir
cp RefSeg1Exact.txt RefSeg2Exact.txt RefSeg3Exact.txt TestInterpolating.testdir/
cp WeaveOut.fits TestInterpolating.testdir/RefWeaveOut.fits
tar zcf TestInterpolating.tar.gz TestInterpolating.testdir/
rm -rf TestInterpolating.testdir/

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
#   octapps_run WeaveFstatMismatch --setup-file=TestInterpolating.testdir/WeaveSetup.fits --spindowns=1 --semi-max-mismatch=5.5 --coh-max-mismatch=0.3 --printarg=meanOfHist
mean_mu=0.47315
awk "BEGIN { print mu = ( ${coh2F_exact} - ${coh2F_loud} ) / ${coh2F_exact}; exit ( mu < ${mean_mu} ? 0 : 1 ) }"
set +x
echo
