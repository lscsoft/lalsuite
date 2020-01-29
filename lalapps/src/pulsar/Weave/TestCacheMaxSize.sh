# Perform an interpolating search without/with a maximum cache size, and check for consistent results

export LAL_FSTAT_FFT_PLAN_MODE=ESTIMATE

for setup in short long; do

    case ${setup} in

        short)
            verb="Perform"
            weave_setup_options="--segment-count=3"
            weave_sft_options="--rand-seed=3456 --sft-timebase=1800 --sft-noise-sqrtSX=1,1 --sft-timestamps-files=timestamps-1.txt,timestamps-2.txt"
            weave_search_options="--alpha=0.9/1.4 --delta=-1.2/2.3 --freq=50.5/0.01 --f1dot=-1.5e-9,0 --semi-max-mismatch=5 --coh-max-mismatch=0.3"
            weave_cache_options="--cache-max-size=2"
            weave_recomp_threshold=0.0
            ;;

        long)
            verb="Simulate"
            weave_setup_options="--segment-count=3 --segment-gap=11130000"
            weave_sft_options=
            weave_search_options="--simulate-search --alpha=2.3/0.9 --delta=-1.2/2.3 --freq=50.5/0.01 --f1dot=-5e-11,0 --semi-max-mismatch=6 --coh-max-mismatch=0.3"
            weave_cache_options="--cache-max-size=25 --cache-all-gc"
            weave_recomp_threshold=0.0
            ;;

        *)
            echo "$0: unknown setup '${setup}'"
            exit 1

    esac

    echo "=== Setup '${setup}': Create search setup with ${weave_setup_options} ==="
    set -x
    lalapps_WeaveSetup --first-segment=1122332211/90000 ${weave_setup_options} --detectors=H1,L1 --output-file=WeaveSetup.fits
    lalapps_fits_overview WeaveSetup.fits
    set +x
    echo

    echo "=== Setup '${setup}': Restrict timestamps to segment list in WeaveSetup.fits ==="
    set -x
    lalapps_fits_table_list 'WeaveSetup.fits[segments][col c1=start_s; col2=end_s]' \
        | awk 'BEGIN { print "/^#/ { print }" } /^#/ { next } { printf "%i <= $1 && $1 <= %i { print }\n", $1, $2 + 1 }' > timestamp-filter.awk
    awk -f timestamp-filter.awk all-timestamps-1.txt > timestamps-1.txt
    awk -f timestamp-filter.awk all-timestamps-2.txt > timestamps-2.txt
    set +x
    echo

    echo "=== Setup '${setup}': ${verb} interpolating search without a maximum cache size ==="
    set -x
    lalapps_Weave --cache-max-size=0 --output-file=WeaveOutNoMax.fits \
        --toplists=all --toplist-limit=2321 --segment-info --setup-file=WeaveSetup.fits \
        ${weave_sft_options} ${weave_search_options}
    lalapps_fits_overview WeaveOutNoMax.fits
    set +x
    echo

    echo "=== Setup '${setup}': Check for non-singular semicoherent dimensions ==="
    set -x
    semi_ntmpl_prev=1
    for dim in SSKYA SSKYB NU1DOT NU0DOT; do
        semi_ntmpl=`lalapps_fits_header_getval "WeaveOutNoMax.fits[0]" "NSEMITMPL ${dim}" | tr '\n\r' '  ' | awk 'NF == 1 {printf "%d", $1}'`
        expr ${semi_ntmpl} '>' ${semi_ntmpl_prev}
        semi_ntmpl_prev=${semi_ntmpl}
    done
    set +x
    echo

    echo "=== Setup '${setup}': ${verb} interpolating search with a maximum cache size ==="
    set -x
    lalapps_Weave ${weave_cache_options} --output-file=WeaveOutMax.fits \
        --toplists=all --toplist-limit=2321 --segment-info --setup-file=WeaveSetup.fits \
        ${weave_sft_options} ${weave_search_options}
    lalapps_fits_overview WeaveOutMax.fits
    set +x
    echo

    echo "=== Setup '${setup}': Check that number of coherent templates are equal ==="
    set -x
    coh_ntmpl_no_max=`lalapps_fits_header_getval "WeaveOutNoMax.fits[0]" 'NCOHTPL' | tr '\n\r' '  ' | awk 'NF == 1 {printf "%d", $1}'`
    coh_ntmpl_max=`lalapps_fits_header_getval "WeaveOutMax.fits[0]" 'NCOHTPL' | tr '\n\r' '  ' | awk 'NF == 1 {printf "%d", $1}'`
    expr ${coh_ntmpl_no_max} '=' ${coh_ntmpl_max}
    set +x
    echo

    echo "=== Setup '${setup}': Check that without a maximum cache number of recomputed results is below tolerance ==="
    set -x
    coh_nres_no_max=`lalapps_fits_header_getval "WeaveOutNoMax.fits[0]" 'NCOHRES' | tr '\n\r' '  ' | awk 'NF == 1 {printf "%d", $1}'`
    awk "BEGIN { print recomp = ( ${coh_nres_no_max} - ${coh_ntmpl_no_max} ) / ${coh_ntmpl_no_max}; exit ( recomp <= ${weave_recomp_threshold} ? 0 : 1 ) }"
    set +x
    echo

    echo "=== Setup '${setup}': Check that with a maximum cache number of recomputed results is above tolerance ==="
    set -x
    coh_nres_max=`lalapps_fits_header_getval "WeaveOutMax.fits[0]" 'NCOHRES' | tr '\n\r' '  ' | awk 'NF == 1 {printf "%d", $1}'`
    awk "BEGIN { print recomp = ( ${coh_nres_max} - ${coh_ntmpl_max} ) / ${coh_ntmpl_max}; exit ( recomp > ${weave_recomp_threshold} ? 0 : 1 ) }"
    set +x
    echo

    case ${setup} in

        short)
            echo "=== Setup '${setup}': Compare F-statistics from lalapps_Weave without/with a maximum cache size ==="
            set -x
            env LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},info" lalapps_WeaveCompare --setup-file=WeaveSetup.fits --result-file-1=WeaveOutNoMax.fits --result-file-2=WeaveOutMax.fits
            set +x
            echo
            ;;

        *)
            ;;

    esac

done
