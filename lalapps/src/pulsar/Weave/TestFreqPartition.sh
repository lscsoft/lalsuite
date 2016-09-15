# Perform a non-interpolating search without/with frequency partitions, and check for consistent results

echo "=== Create search setup with 3 segments ==="
set -x
${builddir}/lalapps_WeaveSetup --first-segment=1122332211/90000 --segment-count=3 --detectors=H1,L1 --output-file=WeaveSetup.fits
set +x
echo

echo "=== Perform non-interpolating search without frequency partitions ==="
set -x
${builddir}/lalapps_Weave --freq-partitions=1 --output-file=WeaveOutNoPart.fits \
    --output-toplist-limit=0 --output-misc-info --setup-file=WeaveSetup.fits \
    --sft-timebase=1800 --sft-noise-psd=1,1 --sft-noise-rand-seed=3456 \
    --sft-timestamps-files=${srcdir}/timestamps-irregular.txt,${srcdir}/timestamps-regular.txt \
    --injections="{refTime=1122332211; h0=0.5; cosi=0.2; psi=0.4; phi0=0.1; Alpha=2.72; Delta=-0.38; Freq=50.5; f1dot=-1e-9}" \
    --alpha=2.72 --delta=-0.38 --freq=50/2e-2 --f1dot=-1e-9 --semi-max-mismatch=0.5 --interpolation=no --Fstat-method=DemodBest
set +x
echo

echo "=== Perform non-interpolating search with frequency partitions ==="
set -x
${builddir}/lalapps_Weave --freq-partitions=11 --output-file=WeaveOutPart.fits \
    --output-toplist-limit=0 --output-misc-info --setup-file=WeaveSetup.fits \
    --sft-timebase=1800 --sft-noise-psd=1,1 --sft-noise-rand-seed=3456 \
    --sft-timestamps-files=${srcdir}/timestamps-irregular.txt,${srcdir}/timestamps-regular.txt \
    --injections="{refTime=1122332211; h0=0.5; cosi=0.2; psi=0.4; phi0=0.1; Alpha=2.72; Delta=-0.38; Freq=50.5; f1dot=-1e-9}" \
    --alpha=2.72 --delta=-0.38 --freq=50/2e-2 --f1dot=-1e-9 --semi-max-mismatch=0.5 --interpolation=no --Fstat-method=DemodBest
set +x
echo

echo "=== Compare F-statistics from lalapps_Weave without/with frequency partitions ==="
set -x
LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},info"
${builddir}/lalapps_WeaveCompare --setup-file=WeaveSetup.fits --output-file-1=WeaveOutNoPart.fits --output-file-2=WeaveOutPart.fits
set +x
echo
