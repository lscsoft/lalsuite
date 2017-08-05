# Perform a fully-coherent search of a single segment, and compare F-statistics to reference results

echo "=== Create single-segment search setup ==="
set -x
${builddir}/lalapps_WeaveSetup --first-segment=1122332211/90000 --detectors=H1,L1 --output-file=WeaveSetup.fits
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
inject_params="Alpha=5.4; Delta=1.1; Freq=55.5; f1dot=-1e-9"
${injdir}/lalapps_Makefakedata_v5 --randSeed=1234 --fmin=55.0 --Band=1.0 \
    --injectionSources="{refTime=${ref_time}; h0=0.5; cosi=0.2; psi=0.4; phi0=0.1; ${inject_params}}" \
    --Tsft=1800 --outSingleSFT --outSFTdir=. --IFOs=H1,L1 --sqrtSX=1,1 \
    --timestampsFiles=timestamps-1.txt,timestamps-2.txt
set +x
echo

echo "=== Perform fully-coherent search ==="
set -x
${builddir}/lalapps_Weave --output-file=WeaveOut.fits \
    --toplists=mean2F --toplist-limit=0 --per-detector --misc-info \
    --setup-file=WeaveSetup.fits --sft-files='*.sft' --Fstat-method=DemodBest \
    --freq=55.5~1e-4 --f1dot=-2e-9,0 --semi-max-mismatch=9
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
${fitsdir}/lalapps_fits_header_getval "WeaveOut.fits[0]" 'NCOHRES' > tmp
coh_nres=`cat tmp | xargs printf "%d"`
${fitsdir}/lalapps_fits_header_getval "WeaveOut.fits[0]" 'NSEMIRES' > tmp
semi_nres=`cat tmp | xargs printf "%d"`
expr ${coh_nres} '=' ${semi_nres}
set +x
echo

echo "=== Check that no results were recomputed ==="
set -x
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[per_seg_info][col coh_nrecomp]" > tmp
coh_nrecomp=`cat tmp | sed "/^#/d" | xargs printf "%d"`
expr ${coh_nrecomp} '=' 0
set +x
echo

### Make updating reference results a little easier ###
mkdir TestSingleSegment.testdir
cp RefExact.txt TestSingleSegment.testdir/
cp WeaveOut.fits TestSingleSegment.testdir/RefWeaveOut.fits
tar zcvf TestSingleSegment.tar.gz TestSingleSegment.testdir/
rm -rf TestSingleSegment.testdir/

echo "=== Compare F-statistics from lalapps_Weave to reference results ==="
set -x
env LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},info" ${builddir}/lalapps_WeaveCompare --setup-file=WeaveSetup.fits --result-file-1=WeaveOut.fits --result-file-2=RefWeaveOut.fits
set +x
echo

echo "=== Compare F-statistic at exact injected signal parameters with loudest F-statistic found by lalapps_Weave ==="
set -x
coh2F_exact=`cat RefExact.txt | sed -n '/^[^%]/p' | awk '{print $7}'`
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean2F_toplist][col c1=mean2F][#row == 1]" > tmp
coh2F_loud=`cat tmp | sed "/^#/d" | xargs printf "%g"`
# Value of 'mean_mu' was calculated by:
#   octapps_run WeaveFstatMismatch --setup-file=TestSingleSegment.testdir/WeaveSetup.fits --spindowns=1 --semi-max-mismatch=9 --coh-max-mismatch=0 --printarg=meanOfHist
mean_mu=0.58305
awk "BEGIN { print mu = ( ${coh2F_exact} - ${coh2F_loud} ) / ${coh2F_exact}; exit ( mu < ${mean_mu} ? 0 : 1 ) }"
set +x
echo
