# Perform an interpolating search without/with a maximum cache size, and check for consistent results

export LAL_FSTAT_FFT_PLAN_MODE=ESTIMATE

for setup in short long; do

    case ${setup} in

        short)
            verb="Perform"
            weave_setup_options="--segment-count=3"
            weave_sft_options="--rand-seed=3456 --sft-timebase=1800 --sft-noise-psd=1,1 --sft-timestamps-files=timestamps-1.txt,timestamps-2.txt"
            weave_search_options="--alpha=0.9/1.4 --delta=-1.2/2.3 --freq=50.5/1e-4 --f1dot=-1.5e-9,0 --semi-max-mismatch=5 --coh-max-mismatch=0.3"
            ;;

        long)
            verb="Simulate"
            weave_setup_options="--segment-count=3 --segment-gap=11130000"
            weave_sft_options=
            weave_search_options="--simulate-search --alpha=0.9/0.05 --delta=-1.2/0.1 --freq=50.5/1e-5 --f1dot=-5e-11,0 --semi-max-mismatch=6 --coh-max-mismatch=0.3"
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

    echo "=== Setup '${setup}': ${verb} interpolating search without a maximum cache size ==="
    set -x
    ${builddir}/lalapps_Weave --cache-max-size=0 --output-file=WeaveOutNoMax.fits \
        --toplists=all --toplist-limit=2321 --misc-info --setup-file=WeaveSetup.fits \
        ${weave_sft_options} ${weave_search_options}
    set +x
    echo

    echo "=== Setup '${setup}': Check for non-singular semicoherent dimensions ==="
    set -x
    semi_ntmpl_prev=1
    for dim in SSKYA SSKYB NU1DOT NU0DOT; do
        ${fitsdir}/lalapps_fits_header_getval "WeaveOutNoMax.fits[0]" "NSEMITMPL ${dim}" > tmp
        semi_ntmpl=`cat tmp | xargs printf "%d"`
        expr ${semi_ntmpl} '>' ${semi_ntmpl_prev}
        semi_ntmpl_prev=${semi_ntmpl}
    done
    set +x
    echo

    echo "=== Setup '${setup}': ${verb} interpolating search with a maximum cache size ==="
    set -x
    ${builddir}/lalapps_Weave --cache-max-size=10 --output-file=WeaveOutMax.fits \
        --toplists=all --toplist-limit=2321 --misc-info --setup-file=WeaveSetup.fits \
        ${weave_sft_options} ${weave_search_options}
    set +x
    echo

    for seg in 1 2 3; do

        echo "=== Setup '${setup}': Segment #${seg}: Check that number of computed coherent results are equal ==="
        set -x
        ${fitsdir}/lalapps_fits_table_list "WeaveOutNoMax.fits[per_seg_info][col coh_n1comp][#row == ${seg}]" > tmp
        coh_n1comp_no_max=`cat tmp | sed "/^#/d" | xargs printf "%d"`
        ${fitsdir}/lalapps_fits_table_list "WeaveOutMax.fits[per_seg_info][col coh_n1comp][#row == ${seg}]" > tmp
        coh_n1comp_max=`cat tmp | sed "/^#/d" | xargs printf "%d"`
        expr ${coh_n1comp_no_max} '=' ${coh_n1comp_max}
        set +x
        echo

        echo "=== Setup '${setup}': Segment #${seg}: Check that search without a maximum cache size did not recompute results ==="
        set -x
        ${fitsdir}/lalapps_fits_table_list "WeaveOutNoMax.fits[per_seg_info][col coh_nrecomp][#row == ${seg}]" > tmp
        coh_nrecomp_no_max=`cat tmp | sed "/^#/d" | xargs printf "%d"`
        expr ${coh_nrecomp_no_max} '=' 0
        set +x
        echo

        echo "=== Setup '${setup}': Segment #${seg}: Check that search with a maximum cache size did recompute results ==="
        set -x
        ${fitsdir}/lalapps_fits_table_list "WeaveOutMax.fits[per_seg_info][col coh_nrecomp][#row == ${seg}]" > tmp
        coh_nrecomp_max=`cat tmp | sed "/^#/d" | xargs printf "%d"`
        expr ${coh_nrecomp_max} '>' 0
        set +x
        echo

    done

    case ${setup} in

        short)
            echo "=== Setup '${setup}': Compare F-statistics from lalapps_Weave without/with a maximum cache size ==="
            set -x
            env LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},info" ${builddir}/lalapps_WeaveCompare --setup-file=WeaveSetup.fits --result-file-1=WeaveOutNoMax.fits --result-file-2=WeaveOutMax.fits
            set +x
            echo
            ;;

        *)
            ;;

    esac

done
