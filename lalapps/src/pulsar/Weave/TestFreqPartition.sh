# Perform an interpolating search without/with frequency partitions, and check for consistent results

echo "=== Create search setup with 3 segments ==="
set -x
${builddir}/lalapps_WeaveSetup --first-segment=1122332211/90000 --segment-count=3 --detectors=H1,L1 --output-file=WeaveSetup.fits
set +x
echo

echo "=== Perform non-interpolating search without frequency partitions ==="
set -x
${builddir}/lalapps_Weave --freq-partitions=1 --output-file=WeaveOutNoPart.fits \
    --output-toplist-limit=5000 --output-misc-info --setup-file=WeaveSetup.fits \
    --sft-timebase=1800 --sft-noise-psd=1,1 --sft-noise-rand-seed=3456 \
    --sft-timestamps-files=${srcdir}/timestamps-irregular.txt,${srcdir}/timestamps-regular.txt \
    --injections="{h0=0.5; cosi=0.2; psi=0.4; phi0=0.1; Alpha=2.72; Delta=-0.38; Freq=50.5; f1dot=-1e-9}" \
    --alpha=2.72/0.05 --delta=-0.38/0.05 --freq=50/1e-4 --f1dot=-1e-8,0 --semi-max-mismatch=0.5 --coh-max-mismatch=0.4
set +x
echo

echo "=== Perform non-interpolating search with frequency partitions ==="
set -x
${builddir}/lalapps_Weave --freq-partitions=5 --output-file=WeaveOutPart.fits \
    --output-toplist-limit=5000 --output-misc-info --setup-file=WeaveSetup.fits \
    --sft-timebase=1800 --sft-noise-psd=1,1 --sft-noise-rand-seed=3456 \
    --sft-timestamps-files=${srcdir}/timestamps-irregular.txt,${srcdir}/timestamps-regular.txt \
    --injections="{h0=0.5; cosi=0.2; psi=0.4; phi0=0.1; Alpha=2.72; Delta=-0.38; Freq=50.5; f1dot=-1e-9}" \
    --alpha=2.72/0.05 --delta=-0.38/0.05 --freq=50/1e-4 --f1dot=-1e-8,0 --semi-max-mismatch=0.5 --coh-max-mismatch=0.4
set +x
echo

for seg in 1 2 3; do

    echo "=== Check that no results were recomputed ==="
    set -x
    ${fitsdir}/lalapps_fits_table_list "WeaveOutNoPart.fits[per_seg_info][col coh_total_recomp][#row == ${seg}]" > tmp
    coh_total_recomp_no_part=`cat tmp | sed "/^#/d" | xargs printf "%d"`
    [ ${coh_total_recomp_no_part} -eq 0 ]
    ${fitsdir}/lalapps_fits_table_list "WeaveOutPart.fits[per_seg_info][col coh_total_recomp][#row == ${seg}]" > tmp
    coh_nrecomp_part=`cat tmp | sed "/^#/d" | xargs printf "%d"`
    [ ${coh_nrecomp_part} -eq 0 ]
    set +x
    echo

done

echo "=== Compare F-statistics from lalapps_Weave without/with frequency partitions ==="
set -x
LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},info"
${builddir}/lalapps_WeaveCompare --setup-file=WeaveSetup.fits --output-file-1=WeaveOutNoPart.fits --output-file-2=WeaveOutPart.fits
set +x
echo
