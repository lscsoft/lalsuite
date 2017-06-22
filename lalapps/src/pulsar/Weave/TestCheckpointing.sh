# Perform an interpolating search without/with checkpointing, and check for consistent results

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

echo "=== Generate SFTs ==="
set -x
${injdir}/lalapps_Makefakedata_v5 --randSeed=3456 --fmin=49.5 --Band=2.0 --Tsft=1800 \
    --outSingleSFT --outSFTdir=. --IFOs=H1,L1 --sqrtSX=1,1 \
    --timestampsFiles=timestamps-1.txt,timestamps-2.txt
set +x
echo

echo "=== Perform interpolating search without checkpointing ==="
set -x
${builddir}/lalapps_Weave --output-file=WeaveOutNoCkpt.fits \
    --toplists=all --toplist-limit=2321 --per-detector --per-segment --setup-file=WeaveSetup.fits --sft-files='*.sft' \
    --sky-patch-count=3 --sky-patch-index=0 --freq=50/1e-4 --f1dot=-1e-9,0 --semi-max-mismatch=5 --coh-max-mismatch=0.4
set +x
echo

echo "=== Check for non-singular semicoherent dimensions ==="
set -x
semi_ntmpl_prev=1
for dim in SSKYA SSKYB NU1DOT NU0DOT; do
    ${fitsdir}/lalapps_fits_header_getval "WeaveOutNoCkpt.fits[0]" "NSEMITMPL ${dim}" > tmp
    semi_ntmpl=`cat tmp | xargs printf "%d"`
    expr ${semi_ntmpl} '>' ${semi_ntmpl_prev}
    semi_ntmpl_prev=${semi_ntmpl}
done
set +x
echo

echo "=== Perform interpolating search with checkpointing ==="
echo "--- Start to first checkpoint ---"
set -x
rm -f WeaveCkpt.fits
${builddir}/lalapps_Weave --output-file=WeaveOutCkpt.fits --ckpt-output-file=WeaveCkpt.fits --ckpt-output-exit=0.22 \
    --toplists=all --toplist-limit=2321 --per-detector --per-segment --setup-file=WeaveSetup.fits --sft-files='*.sft' \
    --sky-patch-count=3 --sky-patch-index=0 --freq=50/1e-4 --f1dot=-1e-9,0 --semi-max-mismatch=5 --coh-max-mismatch=0.4
set +x
echo "--- First to second checkpoint ---"
set -x
${builddir}/lalapps_Weave --output-file=WeaveOutCkpt.fits --ckpt-output-file=WeaveCkpt.fits --ckpt-output-exit=0.63 \
    --toplists=all --toplist-limit=2321 --per-detector --per-segment --setup-file=WeaveSetup.fits --sft-files='*.sft' \
    --sky-patch-count=3 --sky-patch-index=0 --freq=50/1e-4 --f1dot=-1e-9,0 --semi-max-mismatch=5 --coh-max-mismatch=0.4
set +x
echo "--- Second checkpoint to end ---"
set -x
${builddir}/lalapps_Weave --output-file=WeaveOutCkpt.fits --ckpt-output-file=WeaveCkpt.fits \
    --toplists=all --toplist-limit=2321 --per-detector --per-segment --setup-file=WeaveSetup.fits --sft-files='*.sft' \
    --sky-patch-count=3 --sky-patch-index=0 --freq=50/1e-4 --f1dot=-1e-9,0 --semi-max-mismatch=5 --coh-max-mismatch=0.4
set +x
echo

echo "=== Check number of times output results have been restored from a checkpoint ==="
set -x
${fitsdir}/lalapps_fits_header_getval "WeaveCkpt.fits[0]" CKPTCNT > tmp
ckpt_count=`cat tmp | xargs printf "%d"`
expr ${ckpt_count} '=' 2
${fitsdir}/lalapps_fits_header_getval "WeaveOutCkpt.fits[0]" NUMCKPT > tmp
num_ckpt=`cat tmp | xargs printf "%d"`
expr ${num_ckpt} '=' ${ckpt_count}
set +x
echo

echo "=== Compare F-statistics from lalapps_Weave without/with checkpointing ==="
set -x
env LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},info" ${builddir}/lalapps_WeaveCompare --setup-file=WeaveSetup.fits --result-file-1=WeaveOutNoCkpt.fits --result-file-2=WeaveOutCkpt.fits
set +x
echo
