# Perform an interpolating search without/with part of the search being shortcut, and check for consistent results

echo "=== Create search setup with 4 segments ==="
set -x
${builddir}/lalapps_WeaveSetup --first-segment=1122332211/90000 --segment-count=4 --detectors=H1,L1 --output-file=WeaveSetup.fits
set +x
echo

echo "=== Perform interpolating search ==="
set -x
${builddir}/lalapps_Weave --output-file=WeaveOutNoShort.fits \
    --toplist-limit=5000 --misc-info --setup-file=WeaveSetup.fits \
    --rand-seed=3456 --sft-timebase=1800 --sft-noise-psd=1,1 \
    --sft-timestamps-files=${srcdir}/timestamps-irregular.txt,${srcdir}/timestamps-regular.txt \
    --alpha=0.9/1.4 --delta=-1.2/2.3 --freq=55.5/1e-4 --f1dot=-1.5e-9,0 --semi-max-mismatch=6 --coh-max-mismatch=0.3
set +x
echo

echo "=== Check average number of semicoherent templates per dimension is more than one"
set -x
for dim in SSKYA SSKYB NU0DOT NU1DOT; do
    ${fitsdir}/lalapps_fits_header_getval "WeaveOutNoShort.fits[0]" "SEMIAVG ${dim}" > tmp
    semi_avg_ntmpl_dim=`cat tmp | xargs printf "%d"`
    expr ${semi_avg_ntmpl_dim} '>' 1
done
set +x
echo

echo "=== Perform interpolating search with all computations shortcut ==="
set -x
${builddir}/lalapps_Weave --shortcut-compute --output-file=WeaveOutShortComp.fits \
    --toplist-limit=5000 --misc-info --setup-file=WeaveSetup.fits \
    --rand-seed=3456 --sft-timebase=1800 --sft-noise-psd=1,1 \
    --sft-timestamps-files=${srcdir}/timestamps-irregular.txt,${srcdir}/timestamps-regular.txt \
    --alpha=0.9/1.4 --delta=-1.2/2.3 --freq=55.5/1e-4 --f1dot=-1.5e-9,0 --semi-max-mismatch=6 --coh-max-mismatch=0.3
set +x
echo

echo "=== Perform interpolating search with main search loop shortcut ==="
set -x
${builddir}/lalapps_Weave --shortcut-search --output-file=WeaveOutShortSearch.fits \
    --toplist-limit=5000 --misc-info --setup-file=WeaveSetup.fits \
    --rand-seed=3456 --sft-timebase=1800 --sft-noise-psd=1,1 \
    --sft-timestamps-files=${srcdir}/timestamps-irregular.txt,${srcdir}/timestamps-regular.txt \
    --alpha=0.9/1.4 --delta=-1.2/2.3 --freq=55.5/1e-4 --f1dot=-1.5e-9,0 --semi-max-mismatch=6 --coh-max-mismatch=0.3
set +x
echo

echo "=== Check approximate/computed number of semicoherent templates are consistent ==="
set -x
${fitsdir}/lalapps_fits_header_getval "WeaveOutNoShort.fits[0]" 'SEMIAPPX' > tmp
semi_ntmpl_no_short=`cat tmp | xargs printf "%d"`
${fitsdir}/lalapps_fits_header_getval "WeaveOutShortComp.fits[0]" 'SEMIAPPX' > tmp
semi_ntmpl_short_comp=`cat tmp | xargs printf "%d"`
${fitsdir}/lalapps_fits_header_getval "WeaveOutShortSearch.fits[0]" 'SEMIAPPX' > tmp
semi_ntmpl_short_search=`cat tmp | xargs printf "%d"`
expr ${semi_ntmpl_no_short} '=' ${semi_ntmpl_short_comp}
expr ${semi_ntmpl_no_short} '=' ${semi_ntmpl_short_search}
${fitsdir}/lalapps_fits_header_getval "WeaveOutNoShort.fits[0]" 'SEMICOMP' > tmp
semi_ncomp_no_short=`cat tmp | xargs printf "%d"`
${fitsdir}/lalapps_fits_header_getval "WeaveOutShortComp.fits[0]" 'SEMICOMP' > tmp
semi_ncomp_short_comp=`cat tmp | xargs printf "%d"`
expr ${semi_ncomp_no_short} '=' ${semi_ncomp_short_comp}
set +x
echo

echo "=== Check peak memory usage is consistent ==="
set -x
${fitsdir}/lalapps_fits_header_getval "WeaveOutNoShort.fits[0]" 'PEAKMEM' > tmp
peak_mem_no_short=`cat tmp | xargs printf "%g"`
${fitsdir}/lalapps_fits_header_getval "WeaveOutShortComp.fits[0]" 'PEAKMEM' > tmp
peak_mem_short_comp=`cat tmp | xargs printf "%g"`
awk "BEGIN { print x = ${peak_mem_short_comp} / ${peak_mem_no_short}; exit ( ( 0.95 < x && x < 1.05 ) ? 0 : 1 ) }"
set +x
echo

echo "=== Check miscellaneous per-segment information is consistent ==="
set -x
${fitsdir}/lalapps_fits_table_list "WeaveOutNoShort.fits[per_seg_info]" > per_seg_info_no_short
${fitsdir}/lalapps_fits_table_list "WeaveOutShortComp.fits[per_seg_info]" > per_seg_info_short_comp
diff per_seg_info_no_short per_seg_info_short_comp
set +x
echo
