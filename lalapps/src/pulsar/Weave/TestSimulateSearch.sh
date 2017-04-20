# Perform an interpolating search without/with simulation, and check for consistent results

for setup in short mid long; do

    case ${setup} in

        short)
            weave_setup_options="--segment-count=4"
            weave_search_options="--alpha=0.9/1.4 --delta=-1.2/2.3 --freq=55.5/1e-4 --f1dot=-1.5e-9,0 --semi-max-mismatch=6 --coh-max-mismatch=0.3"
            ;;

        mid)
            weave_setup_options="--segment-count=2 --segment-gap=11130000"
            weave_search_options="--alpha=0.9/0.6 --delta=-1.2/1.5 --freq=55.5/1e-4 --f1dot=-1.5e-9,0 --semi-max-mismatch=7 --coh-max-mismatch=0.4"
            ;;

        long)
            weave_setup_options="--segment-count=2 --segment-gap=22350000"
            weave_search_options="--alpha=0.9/0.4 --delta=-1.2/1.3 --freq=55.5/1e-4 --f1dot=-1.5e-9,0 --semi-max-mismatch=5 --coh-max-mismatch=0.2"
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

    echo "=== Setup '${setup}': Perform interpolating search ==="
    set -x
    ${builddir}/lalapps_Weave --output-file=WeaveOutNoSim.fits \
        --toplists=all --toplist-limit=2321 --misc-info --setup-file=WeaveSetup.fits \
        --rand-seed=3456 --sft-timebase=1800 --sft-noise-psd=1,1 \
        --sft-timestamps-files=timestamps-1.txt,timestamps-2.txt \
        ${weave_search_options}
    set +x
    echo

    echo "=== Setup '${setup}': Check average number of semicoherent templates per dimension ==="
    set -x
    for dim in SSKYA SSKYB NU0DOT NU1DOT; do
        ${fitsdir}/lalapps_fits_header_getval "WeaveOutNoSim.fits[0]" "SEMIAVG ${dim}" > tmp
        semi_avg_dim=`cat tmp | xargs printf "%g"`
        awk "BEGIN { exit ( ${semi_avg_dim} > 1 ? 0 : 1 ) }"
    done
    set +x
    echo

    echo "=== Setup '${setup}': Simulate interpolating search with full memory allocation ==="
    set -x
    ${builddir}/lalapps_Weave --simulate-search --output-file=WeaveOutSimFull.fits \
        --toplists=all --toplist-limit=2321 --misc-info --setup-file=WeaveSetup.fits \
        --rand-seed=3456 --sft-timebase=1800 --sft-noise-psd=1,1 \
        --sft-timestamps-files=timestamps-1.txt,timestamps-2.txt \
        ${weave_search_options}
    set +x
    echo

    echo "=== Setup '${setup}': Simulate interpolating search with minimal memory allocation ==="
    set -x
    ${builddir}/lalapps_Weave --simulate-search --output-file=WeaveOutSimMin.fits \
        --toplists=all --toplist-limit=2321 --misc-info --setup-file=WeaveSetup.fits \
        ${weave_search_options}
    set +x
    echo

    echo "=== Setup '${setup}': Check number of coherent frequency blocks ==="
    set -x
    ${fitsdir}/lalapps_fits_header_getval "WeaveOutNoSim.fits[0]" 'NCOHFBK' > tmp
    coh_nfbk_no_sim=`cat tmp | xargs printf "%d"`
    ${fitsdir}/lalapps_fits_header_getval "WeaveOutSimFull.fits[0]" 'NCOHFBK' > tmp
    coh_nfbk_sim_full=`cat tmp | xargs printf "%d"`
    ${fitsdir}/lalapps_fits_header_getval "WeaveOutSimMin.fits[0]" 'NCOHFBK' > tmp
    coh_nfbk_sim_min=`cat tmp | xargs printf "%d"`
    expr ${coh_nfbk_no_sim} '=' ${coh_nfbk_sim_full}
    expr ${coh_nfbk_no_sim} '=' ${coh_nfbk_sim_min}
    set +x
    echo

    echo "=== Setup '${setup}': Check number of coherent results ==="
    set -x
    ${fitsdir}/lalapps_fits_header_getval "WeaveOutNoSim.fits[0]" 'NCOHRES' > tmp
    coh_nres_no_sim=`cat tmp | xargs printf "%d"`
    ${fitsdir}/lalapps_fits_header_getval "WeaveOutSimFull.fits[0]" 'NCOHRES' > tmp
    coh_nres_sim_full=`cat tmp | xargs printf "%d"`
    ${fitsdir}/lalapps_fits_header_getval "WeaveOutSimMin.fits[0]" 'NCOHRES' > tmp
    coh_nres_sim_min=`cat tmp | xargs printf "%d"`
    expr ${coh_nres_no_sim} '=' ${coh_nres_sim_full}
    expr ${coh_nres_no_sim} '=' ${coh_nres_sim_min}
    set +x
    echo

    echo "=== Setup '${setup}': Check number of semicoherent results ==="
    set -x
    ${fitsdir}/lalapps_fits_header_getval "WeaveOutNoSim.fits[0]" 'NSEMIRES' > tmp
    semi_nres_no_sim=`cat tmp | xargs printf "%d"`
    ${fitsdir}/lalapps_fits_header_getval "WeaveOutSimFull.fits[0]" 'NSEMIRES' > tmp
    semi_nres_sim_full=`cat tmp | xargs printf "%d"`
    ${fitsdir}/lalapps_fits_header_getval "WeaveOutSimMin.fits[0]" 'NSEMIRES' > tmp
    semi_nres_sim_min=`cat tmp | xargs printf "%d"`
    expr ${semi_nres_no_sim} '=' ${semi_nres_sim_full}
    expr ${semi_nres_no_sim} '=' ${semi_nres_sim_min}
    set +x
    echo

    echo "=== Setup '${setup}': Check peak memory usage ==="
    set -x
    ${fitsdir}/lalapps_fits_header_getval "WeaveOutNoSim.fits[0]" 'PEAKMEM' > tmp
    peak_mem_no_sim=`cat tmp | xargs printf "%g"`
    ${fitsdir}/lalapps_fits_header_getval "WeaveOutSimFull.fits[0]" 'PEAKMEM' > tmp
    peak_mem_sim_full=`cat tmp | xargs printf "%g"`
    awk "BEGIN { print x = ${peak_mem_sim_full} / ${peak_mem_no_sim}; exit ( ( 0.9 < x && x < 1.0 ) ? 0 : 1 ) }"
    set +x
    echo

    echo "=== Setup '${setup}': Check miscellaneous per-segment information ==="
    set -x
    ${fitsdir}/lalapps_fits_table_list "WeaveOutNoSim.fits[per_seg_info]" > per_seg_info_no_sim
    ${fitsdir}/lalapps_fits_table_list "WeaveOutSimFull.fits[per_seg_info]" > per_seg_info_sim_full
    diff per_seg_info_no_sim per_seg_info_sim_full
    set +x
    echo

done
