# Perform an interpolating search without/with a maximum cache size, and check for consistent results

echo "=== Create search setup with 3 segments ==="
set -x
${builddir}/lalapps_WeaveSetup --first-segment=1122332211/90000 --segment-count=3 --detectors=H1,L1 --output-file=WeaveSetup.fits
set +x
echo

echo "=== Perform interpolating search without a maximum cache size ==="
set -x
${builddir}/lalapps_Weave --cache-max-size=0 \
    --output-file=WeaveOutNoMax.fits --output-max-size=5000 --output-info-per-seg --setup-file=WeaveSetup.fits \
    --sft-detectors=H1,L1 --sft-noise-psd=1,1 --sft-noise-rand-seed=3456 --sft-timebase=1800 \
    --sft-timestamps-files=${srcdir}/timestamps-irregular.txt,${srcdir}/timestamps-regular.txt \
    --injections="{refTime=1122332211; h0=0.5; cosi=0.2; psi=0.4; phi0=0.1; Alpha=2.72; Delta=-0.38; Freq=50.5; f1dot=-1e-9}" \
    --alpha=2.7/0.05 --delta=-0.4/0.05 --freq=50.5/1e-6 --f1dot=-1e-8,0 --semi-max-mismatch=0.6 --coh-max-mismatch=0.3 --Fstat-method=DemodBest
set +x
echo

echo "=== Perform interpolating search with a maximum cache size ==="
set -x
${builddir}/lalapps_Weave --cache-max-size=50 \
    --output-file=WeaveOutMax.fits --output-max-size=5000 --output-info-per-seg --setup-file=WeaveSetup.fits \
    --sft-detectors=H1,L1 --sft-noise-psd=1,1 --sft-noise-rand-seed=3456 --sft-timebase=1800 \
    --sft-timestamps-files=${srcdir}/timestamps-irregular.txt,${srcdir}/timestamps-regular.txt \
    --injections="{refTime=1122332211; h0=0.5; cosi=0.2; psi=0.4; phi0=0.1; Alpha=2.72; Delta=-0.38; Freq=50.5; f1dot=-1e-9}" \
    --alpha=2.7/0.05 --delta=-0.4/0.05 --freq=50.5/1e-6 --f1dot=-1e-8,0 --semi-max-mismatch=0.6 --coh-max-mismatch=0.3 --Fstat-method=DemodBest
set +x
echo

echo "=== Check that coherent template counts are equal, and that WeaveOut{NoMax|Max}.fits {did not|did} recompute results ==="
set -x
for seg in 1 2 3; do
    ${fitsdir}/lalapps_fits_table_list "WeaveOutNoMax.fits[per_seg_info][col coh_total][#row == ${seg}]" > tmp
    coh_total_no_max=`cat tmp | sed "/^#/d" | xargs printf "%d"`
    ${fitsdir}/lalapps_fits_table_list "WeaveOutNoMax.fits[per_seg_info][col coh_total_recomp][#row == ${seg}]" > tmp
    coh_total_recomp_no_max=`cat tmp | sed "/^#/d" | xargs printf "%d"`
    ${fitsdir}/lalapps_fits_table_list "WeaveOutMax.fits[per_seg_info][col coh_total][#row == ${seg}]" > tmp
    coh_total_max=`cat tmp | sed "/^#/d" | xargs printf "%d"`
    ${fitsdir}/lalapps_fits_table_list "WeaveOutMax.fits[per_seg_info][col coh_total_recomp][#row == ${seg}]" > tmp
    coh_total_recomp_max=`cat tmp | sed "/^#/d" | xargs printf "%d"`
    [ ${coh_total_no_max} -eq ${coh_total_no_max} ]
    [ ${coh_total_recomp_no_max} -eq 0 ]
    [ ${coh_total_recomp_max} -gt 0 ]
done
set +x
echo

echo "=== Extract F-statistics from WeaveOut{NoMax|Max}.fits as ASCII table ==="
set -x
${fitsdir}/lalapps_fits_table_list "WeaveOutNoMax.fits[toplist_mean_twoF][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0; c7=DEFNULL(mean_twoF,-999)]" > WeaveFstatsNoMax.txt
${fitsdir}/lalapps_fits_table_list "WeaveOutMax.fits[toplist_mean_twoF][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0; c7=DEFNULL(mean_twoF,-999)]" > WeaveFstatsMax.txt
set +x
echo

echo "=== Compare F-statistics from lalapps_Weave without/with frequency partitions ==="
set -x
${fstatdir}/lalapps_compareFstats --Fname1=WeaveFstatsNoMax.txt --Fname2=WeaveFstatsMax.txt
set +x
echo
