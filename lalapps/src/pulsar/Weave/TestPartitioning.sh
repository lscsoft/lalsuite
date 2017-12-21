# Perform an interpolating search without/with frequency/spindown partitions, and check for consistent results

export LAL_FSTAT_FFT_PLAN_MODE=ESTIMATE

for setup in onefreq short long; do

    case ${setup} in

        onefreq)
            weave_setup_options="--segment-count=3"
            weave_search_options="--alpha=0.2 --delta=-0.4 --freq=47.7 --f1dot=-1.1e-9,0 --semi-max-mismatch=5 --coh-max-mismatch=0.3"
            weave_part_options="--freq-partitions=5"
            weave_recomp_threshold_no_part=0.0
            weave_recomp_threshold_part=0.0
            ;;

        short)
            weave_setup_options="--segment-count=3"
            weave_search_options="--alpha=1.9/1.4 --delta=-1.2/2.3 --freq=49.5/0.01 --f1dot=-1e-9,0 --semi-max-mismatch=6 --coh-max-mismatch=0.3"
            weave_part_options="--freq-partitions=5 --f1dot-partitions=2"
            weave_recomp_threshold_no_part=0.0
            weave_recomp_threshold_part=0.02
            ;;

        long)
            weave_setup_options="--segment-count=3 --segment-gap=11130000"
            weave_search_options="--alpha=0.1/0.5 --delta=-0.2/0.4 --freq=41.5/0.01 --f1dot=-3e-11,0 --semi-max-mismatch=12 --coh-max-mismatch=0.6"
            weave_part_options="--freq-partitions=5"
            weave_recomp_threshold_no_part=0.0
            weave_recomp_threshold_part=0.0
            ;;

        *)
            echo "$0: unknown setup '${setup}'"
            exit 1

    esac

    echo "=== Setup '${setup}': Create search setup with ${weave_setup_options} ==="
    set -x
    ${builddir}/lalapps_WeaveSetup --first-segment=1122332211/90000 ${weave_setup_options} --detectors=H1,L1 --output-file=WeaveSetup.fits
    set +x
    echo

    echo "=== Setup '${setup}': Restrict timestamps to segment list in WeaveSetup.fits ==="
    set -x
    ${fitsdir}/lalapps_fits_table_list 'WeaveSetup.fits[segments][col c1=start_s; col2=end_s]' \
        | awk 'BEGIN { print "/^#/ { print }" } /^#/ { next } { printf "%i <= $1 && $1 <= %i { print }\n", $1, $2 + 1 }' > timestamp-filter.awk
    awk -f timestamp-filter.awk ${srcdir}/timestamps-1.txt > timestamps-1.txt
    awk -f timestamp-filter.awk ${srcdir}/timestamps-2.txt > timestamps-2.txt
    set +x
    echo

    echo "=== Setup '${setup}': Perform interpolating search without frequency/spindown partitions ==="
    set -x
    ${builddir}/lalapps_Weave --output-file=WeaveOutNoPart.fits \
        --toplists=mean2F --toplist-limit=2321 --segment-info --setup-file=WeaveSetup.fits \
        --Fstat-run-med-window=50 \
        --rand-seed=3456 --sft-timebase=1800 --sft-noise-psd=1,1 \
        --sft-timestamps-files=timestamps-1.txt,timestamps-2.txt \
        ${weave_search_options}
    set +x
    echo

    case ${setup} in

        onefreq)
            ;;

        *)
            echo "=== Setup '${setup}': Check for non-singular semicoherent dimensions ==="
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

    esac

    echo "=== Setup '${setup}': Perform interpolating search with frequency/spindown partitions ==="
    set -x
    ${builddir}/lalapps_Weave ${weave_part_options} --output-file=WeaveOutPart.fits \
        --toplists=mean2F --toplist-limit=2321 --segment-info --setup-file=WeaveSetup.fits \
        --Fstat-run-med-window=50 \
        --rand-seed=3456 --sft-timebase=1800 --sft-noise-psd=1,1 \
        --sft-timestamps-files=timestamps-1.txt,timestamps-2.txt \
        ${weave_search_options}
    set +x
    echo

    echo "=== Setup '${setup}': Check that number of recomputed results is below tolerance ==="
    set -x
    ${fitsdir}/lalapps_fits_header_getval "WeaveOutNoPart.fits[0]" 'NCOHRES' > tmp
    coh_nres_no_part=`cat tmp | xargs printf "%d"`
    ${fitsdir}/lalapps_fits_header_getval "WeaveOutNoPart.fits[0]" 'NCOHTPL' > tmp
    coh_ntmpl_no_part=`cat tmp | xargs printf "%d"`
    awk "BEGIN { print recomp = ( ${coh_nres_no_part} - ${coh_ntmpl_no_part} ) / ${coh_ntmpl_no_part}; exit ( recomp <= ${weave_recomp_threshold_no_part} ? 0 : 1 ) }"
    ${fitsdir}/lalapps_fits_header_getval "WeaveOutPart.fits[0]" 'NCOHRES' > tmp
    coh_nres_part=`cat tmp | xargs printf "%d"`
    ${fitsdir}/lalapps_fits_header_getval "WeaveOutPart.fits[0]" 'NCOHTPL' > tmp
    coh_ntmpl_part=`cat tmp | xargs printf "%d"`
    awk "BEGIN { print recomp = ( ${coh_nres_part} - ${coh_ntmpl_part} ) / ${coh_ntmpl_part}; exit ( recomp <= ${weave_recomp_threshold_part} ? 0 : 1 ) }"
    set +x
    echo

    echo "=== Setup '${setup}': Compare F-statistics from lalapps_Weave without/with frequency/spindown partitions ==="
    set -x
    env LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},info" ${builddir}/lalapps_WeaveCompare --setup-file=WeaveSetup.fits --result-file-1=WeaveOutNoPart.fits --result-file-2=WeaveOutPart.fits
    set +x
    echo

done
