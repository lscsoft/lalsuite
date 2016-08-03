# Perform a non-interpolating search without/with frequency partitions, and check for consistent results

echo "=== Create search setup with 3 segments ==="
set -x
${builddir}/lalapps_WeaveSetup --first-segment=1122332211/90000 --segment-count=3 --detectors=H1,L1 --output-file=WeaveSetup.fits
set +x
echo

echo "=== Perform non-interpolating search without frequency partitions ==="
set -x
${builddir}/lalapps_Weave --freq-partitions=1 \
    --output-file=WeaveOutNoPart.fits --output-max-size=0 --output-details --setup-file=WeaveSetup.fits \
    --sft-detectors=H1,L1 --sft-noise-psd=1,1 --sft-noise-rand-seed=3456 --sft-timebase=1800 \
    --sft-timestamps-files=${srcdir}/timestamps-irregular.txt,${srcdir}/timestamps-regular.txt \
    --injections="{refTime=1122332211; h0=0.5; cosi=0.2; psi=0.4; phi0=0.1; Alpha=2.72; Delta=-0.38; Freq=50.5; f1dot=-1e-9}" \
    --alpha=2.72 --delta=-0.38 --freq=50/2e-2 --f1dot=-1e-9 --semi-max-mismatch=0.5 --interpolation=no --Fstat-method=DemodBest
set +x
echo

echo "=== Perform non-interpolating search with frequency partitions ==="
set -x
${builddir}/lalapps_Weave --freq-partitions=11 \
    --output-file=WeaveOutPart.fits --output-max-size=0 --output-details --setup-file=WeaveSetup.fits \
    --sft-detectors=H1,L1 --sft-noise-psd=1,1 --sft-noise-rand-seed=3456 --sft-timebase=1800 \
    --sft-timestamps-files=${srcdir}/timestamps-irregular.txt,${srcdir}/timestamps-regular.txt \
    --injections="{refTime=1122332211; h0=0.5; cosi=0.2; psi=0.4; phi0=0.1; Alpha=2.72; Delta=-0.38; Freq=50.5; f1dot=-1e-9}" \
    --alpha=2.72 --delta=-0.38 --freq=50/2e-2 --f1dot=-1e-9 --semi-max-mismatch=0.5 --interpolation=no --Fstat-method=DemodBest
set +x
echo

echo "=== Extract F-statistics from WeaveOut{NoPart|Part}.fits as ASCII table ==="
set -x
${fitsdir}/lalapps_fits_table_list "WeaveOutNoPart.fits[toplist_mean_twoF][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0; c7=DEFNULL(mean_twoF,-999)]" > WeaveFstatsNoPart.txt
${fitsdir}/lalapps_fits_table_list "WeaveOutPart.fits[toplist_mean_twoF][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0; c7=DEFNULL(mean_twoF,-999)]" > WeaveFstatsPart.txt
set +x
echo

echo "=== Compare F-statistics from lalapps_Weave without/with frequency partitions ==="
set -x
${fstatdir}/lalapps_compareFstats --Fname1=WeaveFstatsNoPart.txt --Fname2=WeaveFstatsPart.txt
set +x
echo
