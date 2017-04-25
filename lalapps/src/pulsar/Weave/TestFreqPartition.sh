# Perform an interpolating search without/with frequency partitions, and check for consistent results

echo "=== Create search setup with 3 segments ==="
set -x
${builddir}/lalapps_WeaveSetup --first-segment=1122332211/90000 --segment-count=3 --detectors=H1,L1 --output-file=WeaveSetup.fits
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

echo "=== Perform non-interpolating search without frequency partitions ==="
set -x
${builddir}/lalapps_Weave --freq-partitions=1 --output-file=WeaveOutNoPart.fits \
    --toplists=all --toplist-limit=2321 --misc-info --setup-file=WeaveSetup.fits \
    --rand-seed=3456 --sft-timebase=1800 --sft-noise-psd=1,1 \
    --sft-timestamps-files=timestamps-1.txt,timestamps-2.txt \
    --alpha=0.9/1.4 --delta=-1.2/2.3 --freq=49.5/1e-3 --f1dot=-1e-9,0 --semi-max-mismatch=6 --coh-max-mismatch=0.3
set +x
echo

echo "=== Check for non-singular semicoherent dimensions ==="
set -x
semi_ntmpl_prev=1
for dim in SSKYA SSKYB NU1DOT NU0DOT; do
    ${fitsdir}/lalapps_fits_header_getval "WeaveOutNoPart.fits[0]" "NSEMITMPL ${dim}" > tmp
    semi_ntmpl=`cat tmp | xargs printf "%d"`
    expr ${semi_ntmpl} '>' ${semi_ntmpl_prev}
    semi_ntmpl_prev=${semi_ntmpl}
done
set +x
echo

echo "=== Perform non-interpolating search with frequency partitions ==="
set -x
${builddir}/lalapps_Weave --freq-partitions=5 --output-file=WeaveOutPart.fits \
    --toplists=all --toplist-limit=2321 --misc-info --setup-file=WeaveSetup.fits \
    --rand-seed=3456 --sft-timebase=1800 --sft-noise-psd=1,1 \
    --sft-timestamps-files=timestamps-1.txt,timestamps-2.txt \
    --alpha=0.9/1.4 --delta=-1.2/2.3 --freq=49.5/1e-3 --f1dot=-1e-9,0 --semi-max-mismatch=6 --coh-max-mismatch=0.3
set +x
echo

for seg in 1 2 3; do

    echo "=== Segment #${seg}: Check that no results were recomputed ==="
    set -x
    ${fitsdir}/lalapps_fits_table_list "WeaveOutNoPart.fits[per_seg_info][col coh_nrecomp][#row == ${seg}]" > tmp
    coh_nrecomp_no_part=`cat tmp | sed "/^#/d" | xargs printf "%d"`
    expr ${coh_nrecomp_no_part} '=' 0
    ${fitsdir}/lalapps_fits_table_list "WeaveOutPart.fits[per_seg_info][col coh_nrecomp][#row == ${seg}]" > tmp
    coh_nrecomp_part=`cat tmp | sed "/^#/d" | xargs printf "%d"`
    expr ${coh_nrecomp_part} '=' 0
    set +x
    echo

done

echo "=== Compare F-statistics from lalapps_Weave without/with frequency partitions ==="
set -x
env LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},info" ${builddir}/lalapps_WeaveCompare --setup-file=WeaveSetup.fits --result-file-1=WeaveOutNoPart.fits --result-file-2=WeaveOutPart.fits
set +x
echo
