if test "${CFITSIO_ENABLED}" = false; then
    echo "Skipping test: requires CFITSIO"
    exit 77
fi

# Perform an interpolating search without/with checkpointing, and check for consistent results

export LAL_FSTAT_FFT_PLAN_MODE=ESTIMATE

echo "=== Create search setup with 3 segments ==="
set -x
lalpulsar_WeaveSetup --first-segment=1122332211/90000 --segment-count=3 --detectors=H1,L1 --output-file=WeaveSetup.fits
lalpulsar_fits_overview WeaveSetup.fits
set +x
echo

echo "=== Restrict timestamps to segment list in WeaveSetup.fits ==="
set -x
lalpulsar_fits_table_list 'WeaveSetup.fits[segments][col c1=start_s; col2=end_s]' \
    | awk 'BEGIN { print "/^#/ { print }" } /^#/ { next } { printf "%i <= $1 && $1 <= %i { print }\n", $1, $2 + 1 }' > timestamp-filter.awk
awk -f timestamp-filter.awk all-timestamps-1.txt > timestamps-1.txt
awk -f timestamp-filter.awk all-timestamps-2.txt > timestamps-2.txt
set +x
echo

echo "=== Generate SFTs ==="
set -x
lalpulsar_Makefakedata_v5 --randSeed=3456 --fmin=49.5 --Band=2.0 --Tsft=1800 \
    --outSingleSFT --outSFTdir=. --IFOs=H1,L1 --sqrtSX=1,1 \
    --timestampsFiles=timestamps-1.txt,timestamps-2.txt
set +x
echo

echo "=== Perform interpolating search without checkpointing ==="
set -x
lalpulsar_Weave --output-file=WeaveOutNoCkpt.fits \
    --toplists=all --toplist-limit=232 --setup-file=WeaveSetup.fits --sft-files='*.sft' --extra-statistics="mean2F_det,sum2F_det,coh2F,coh2F_det" \
    --sky-patch-count=4 --sky-patch-index=0 --freq=50/0.01 --f1dot=-1e-9,0 --semi-max-mismatch=5 --coh-max-mismatch=0.4 --recalc-statistics=all
lalpulsar_fits_overview WeaveOutNoCkpt.fits
set +x
echo

echo "=== Check for non-singular semicoherent dimensions ==="
set -x
semi_ntmpl_prev=1
for dim in SSKYA SSKYB NU1DOT NU0DOT; do
    semi_ntmpl=`lalpulsar_fits_header_getval "WeaveOutNoCkpt.fits[0]" "NSEMITMPL ${dim}" | tr '\n\r' '  ' | awk 'NF == 1 {printf "%d", $1}'`
    expr ${semi_ntmpl} '>' ${semi_ntmpl_prev}
    semi_ntmpl_prev=${semi_ntmpl}
done
set +x
echo

for opt in nopart freqpart f1dotpart allpart; do

    case ${opt} in

        nopart)
            weave_part_options=""
            ;;

        freqpart)
            weave_part_options="--freq-partitions=2"
            ;;

        f1dotpart)
            weave_part_options="--f1dot-partitions=2"
            ;;

        allpart)
            weave_part_options="--freq-partitions=2 --f1dot-partitions=2"
            ;;

        *)
            echo "$0: unknown options '${opt}'"
            exit 1

    esac

    echo "=== Options '${opt}': Perform interpolating search with checkpointing ==="
    echo "--- Start to first checkpoint ---"
    set -x
    rm -f WeaveCkpt.fits
    lalpulsar_Weave ${weave_part_options} --output-file=WeaveOutCkpt.fits --ckpt-output-file=WeaveCkpt.fits --ckpt-output-exit=0.22 \
        --toplists=all --toplist-limit=232 --extra-statistics="mean2F_det,sum2F_det,coh2F,coh2F_det" --setup-file=WeaveSetup.fits --sft-files='*.sft' \
        --sky-patch-count=4 --sky-patch-index=0 --freq=50/0.01 --f1dot=-1e-9,0 --semi-max-mismatch=5 --coh-max-mismatch=0.4 --recalc-statistics=all
    lalpulsar_fits_overview WeaveCkpt.fits
    set +x
    echo "--- First to second checkpoint ---"
    set -x
    lalpulsar_Weave ${weave_part_options} --output-file=WeaveOutCkpt.fits --ckpt-output-file=WeaveCkpt.fits --ckpt-output-exit=0.63 \
        --toplists=all --toplist-limit=232 --extra-statistics="mean2F_det,sum2F_det,coh2F,coh2F_det" --setup-file=WeaveSetup.fits --sft-files='*.sft' \
        --sky-patch-count=4 --sky-patch-index=0 --freq=50/0.01 --f1dot=-1e-9,0 --semi-max-mismatch=5 --coh-max-mismatch=0.4 --recalc-statistics=all
    lalpulsar_fits_overview WeaveCkpt.fits
    set +x
    echo "--- Second checkpoint to end ---"
    set -x
    lalpulsar_Weave ${weave_part_options} --output-file=WeaveOutCkpt.fits --ckpt-output-file=WeaveCkpt.fits \
        --toplists=all --toplist-limit=232 --extra-statistics="mean2F_det,sum2F_det,coh2F,coh2F_det" --setup-file=WeaveSetup.fits --sft-files='*.sft' \
        --sky-patch-count=4 --sky-patch-index=0 --freq=50/0.01 --f1dot=-1e-9,0 --semi-max-mismatch=5 --coh-max-mismatch=0.4 --recalc-statistics=all
    lalpulsar_fits_overview WeaveOutCkpt.fits
    set +x
    echo

    echo "=== Options '${opt}': Check number of times output results have been restored from a checkpoint ==="
    set -x
    ckpt_count=`lalpulsar_fits_header_getval "WeaveCkpt.fits[0]" CKPTCNT | tr '\n\r' '  ' | awk 'NF == 1 {printf "%d", $1}'`
    expr ${ckpt_count} '=' 2
    num_ckpt=`lalpulsar_fits_header_getval "WeaveOutCkpt.fits[0]" NUMCKPT | tr '\n\r' '  ' | awk 'NF == 1 {printf "%d", $1}'`
    expr ${num_ckpt} '=' ${ckpt_count}
    set +x
    echo

    echo "=== Options '${opt}': Compare F-statistics from lalpulsar_Weave without/with checkpointing ==="
    set -x
    lalpulsar_WeaveCompare --setup-file=WeaveSetup.fits --result-file-1=WeaveOutNoCkpt.fits --result-file-2=WeaveOutCkpt.fits --sort-by-semi-phys
    set +x
    echo

done
