# Perform an interpolating search without/with checkpointing, and check for consistent results

echo "=== Generate SFTs ==="
set -x
${injdir}/lalapps_Makefakedata_v5 --randSeed=3456 --fmin=49.5 --Band=2.0 --Tsft=1800 \
    --injectionSources="{refTime=1122332211; h0=0.5; cosi=0.2; psi=0.4; phi0=0.1; Alpha=2.72; Delta=-0.38; Freq=50.5; f1dot=-1e-9}" \
    --outSingleSFT --outSFTdir=. --IFOs=H1,L1 --sqrtSX=1,1 \
    --timestampsFiles=${srcdir}/timestamps-irregular.txt,${srcdir}/timestamps-regular.txt
set +x
echo

echo "=== Create search setup with 3 segments ==="
set -x
${builddir}/lalapps_WeaveSetup --first-segment=1122332211/90000 --segment-count=3 --detectors=H1,L1 --output-file=WeaveSetup.fits
set +x
echo

echo "=== Perform interpolating search without checkpointing ==="
set -x
${builddir}/lalapps_Weave --output-file=WeaveOutNoCkpt.fits \
    --output-toplist-limit=5000 --output-per-detector --output-per-segment --setup-file=WeaveSetup.fits --sft-files='*.sft' \
    --alpha=2.72/0.05 --delta=-0.38/0.05 --freq=50/1e-4 --f1dot=-1e-8,0 --semi-max-mismatch=0.5 --coh-max-mismatch=0.4 --Fstat-method=DemodBest
set +x
echo

echo "=== Perform interpolating search with checkpointing ==="
echo "--- Start to first checkpoint ---"
set -x
rm -f WeaveCkpt.fits
${builddir}/lalapps_Weave --output-file=WeaveOutCkpt.fits --ckpt-output-file=WeaveCkpt.fits --ckpt-output-pc-exit=22 \
    --output-toplist-limit=5000 --output-per-detector --output-per-segment --setup-file=WeaveSetup.fits --sft-files='*.sft' \
    --alpha=2.72/0.05 --delta=-0.38/0.05 --freq=50/1e-4 --f1dot=-1e-8,0 --semi-max-mismatch=0.5 --coh-max-mismatch=0.4 --Fstat-method=DemodBest
set +x
echo "--- First to second checkpoint ---"
set -x
${builddir}/lalapps_Weave --output-file=WeaveOutCkpt.fits --ckpt-output-file=WeaveCkpt.fits --ckpt-output-pc=63 \
    --output-toplist-limit=5000 --output-per-detector --output-per-segment --setup-file=WeaveSetup.fits --sft-files='*.sft' \
    --alpha=2.72/0.05 --delta=-0.38/0.05 --freq=50/1e-4 --f1dot=-1e-8,0 --semi-max-mismatch=0.5 --coh-max-mismatch=0.4 --Fstat-method=DemodBest
set +x
echo "--- Second checkpoint to end ---"
set -x
${builddir}/lalapps_Weave --output-file=WeaveOutCkpt.fits --ckpt-output-file=WeaveCkpt.fits \
    --output-toplist-limit=5000 --output-per-detector --output-per-segment --setup-file=WeaveSetup.fits --sft-files='*.sft' \
    --alpha=2.72/0.05 --delta=-0.38/0.05 --freq=50/1e-4 --f1dot=-1e-8,0 --semi-max-mismatch=0.5 --coh-max-mismatch=0.4 --Fstat-method=DemodBest
set +x
echo

echo "=== Check number of times output results have been restored from a checkpoint ==="
set -x
${fitsdir}/lalapps_fits_header_list "WeaveCkpt.fits[0]" > tmp
ckpt_count=`cat tmp | sed -n "s|/.*$||;s|^CKPTCNT = ||p" | xargs printf "%d"`
[ ${ckpt_count} -eq 2 ]
${fitsdir}/lalapps_fits_header_list "WeaveOutCkpt.fits[0]" > tmp
num_ckpt=`cat tmp | sed -n "s|/.*$||;s|^NUMCKPT = ||p" | xargs printf "%d"`
[ ${num_ckpt} -eq ${ckpt_count} ]
set +x
echo

echo "=== Compare F-statistics from lalapps_Weave without/with checkpointing ==="
set -x
LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},info"
${builddir}/lalapps_WeaveCompare --setup-file=WeaveSetup.fits --output-file-1=WeaveOutNoCkpt.fits --output-file-2=WeaveOutCkpt.fits
set +x
echo
