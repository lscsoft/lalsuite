# Perform an interpolating search without/with a maximum cache size, and check for consistent results

echo "=== Create search setup with 4 segments ==="
set -x
${builddir}/lalapps_WeaveSetup --first-segment=1122332211/90000 --segment-count=4 --detectors=H1,L1 --output-file=WeaveSetup.fits
set +x
echo

echo "=== Perform interpolating search ==="
set -x
${builddir}/lalapps_Weave --output-file=WeaveOutNoSim.fits \
    --output-toplist-limit=5000 --output-misc-info --setup-file=WeaveSetup.fits \
    --sft-timebase=1800 --sft-noise-psd=1,1 --sft-noise-rand-seed=3456 \
    --sft-timestamps-files=${srcdir}/timestamps-irregular.txt,${srcdir}/timestamps-regular.txt \
    --injections="{refTime=1122332211; h0=0.5; cosi=0.2; psi=0.4; phi0=0.1; Alpha=2.72; Delta=-0.38; Freq=50.5; f1dot=-1e-9}" \
    --alpha=2.7/0.05 --delta=-0.4/0.05 --freq=50.5/1e-6 --f1dot=-1e-8,0 --semi-max-mismatch=0.6 --coh-max-mismatch=0.3 --Fstat-method=DemodBest
set +x
echo

echo "=== Simulate interpolating search with full memory usage ==="
set -x
${builddir}/lalapps_Weave --simulate-search --output-file=WeaveOutSimFull.fits \
    --output-toplist-limit=5000 --output-misc-info --setup-file=WeaveSetup.fits \
    --sft-timebase=1800 --sft-noise-psd=1,1 --sft-noise-rand-seed=3456 \
    --sft-timestamps-files=${srcdir}/timestamps-irregular.txt,${srcdir}/timestamps-regular.txt \
    --injections="{refTime=1122332211; h0=0.5; cosi=0.2; psi=0.4; phi0=0.1; Alpha=2.72; Delta=-0.38; Freq=50.5; f1dot=-1e-9}" \
    --alpha=2.7/0.05 --delta=-0.4/0.05 --freq=50.5/1e-6 --f1dot=-1e-8,0 --semi-max-mismatch=0.6 --coh-max-mismatch=0.3 --Fstat-method=DemodBest
set +x
echo

echo "=== Simulate interpolating search with minimal memory usage ==="
set -x
${builddir}/lalapps_Weave --simulate-search --output-file=WeaveOutSimMin.fits \
    --output-toplist-limit=5000 --output-misc-info --setup-file=WeaveSetup.fits \
    --alpha=2.7/0.05 --delta=-0.4/0.05 --freq=50.5/1e-6 --f1dot=-1e-8,0 --semi-max-mismatch=0.6 --coh-max-mismatch=0.3 --Fstat-method=DemodBest
set +x
echo

echo "=== Check number of semicoherent templates are the same ==="
set -x
${fitsdir}/lalapps_fits_header_list "WeaveOutNoSim.fits[0]" > tmp
semi_tot_no_sim=`cat tmp | sed -n "s|/.*$||;s|^SEMITOT = ||p" | xargs printf "%d"`
${fitsdir}/lalapps_fits_header_list "WeaveOutSimFull.fits[0]" > tmp
semi_tot_sim_full=`cat tmp | sed -n "s|/.*$||;s|^SEMITOT = ||p" | xargs printf "%d"`
${fitsdir}/lalapps_fits_header_list "WeaveOutSimMin.fits[0]" > tmp
semi_tot_sim_min=`cat tmp | sed -n "s|/.*$||;s|^SEMITOT = ||p" | xargs printf "%d"`
[ ${semi_tot_no_sim} -eq ${semi_tot_sim_full} ]
[ ${semi_tot_no_sim} -eq ${semi_tot_sim_min} ]
set +x
echo

echo "=== Check peak memory usage is consistent ==="
set -x
${fitsdir}/lalapps_fits_header_list "WeaveOutNoSim.fits[0]" > tmp
peak_mem_no_sim=`cat tmp | sed -n "s|/.*$||;s|^PEAKMEM = ||p" | xargs printf "%.5g"`
${fitsdir}/lalapps_fits_header_list "WeaveOutSimFull.fits[0]" > tmp
peak_mem_sim_full=`cat tmp | sed -n "s|/.*$||;s|^PEAKMEM = ||p" | xargs printf "%.5g"`
${fitsdir}/lalapps_fits_header_list "WeaveOutSimMin.fits[0]" > tmp
peak_mem_sim_min=`cat tmp | sed -n "s|/.*$||;s|^PEAKMEM = ||p" | xargs printf "%.5g"`
[ `echo ${peak_mem_no_sim} ${peak_mem_sim_full} | awk '{ print 0.95 < $2/$1 && $2/$1 < 1.05 }'` -eq 1 ]
[ `echo ${peak_mem_no_sim} ${peak_mem_sim_min} | awk '{ print $2/$1 < 0.6 }'` -eq 1 ]
set +x
echo

echo "=== Check miscellaneous per-segment information is consistent ==="
set -x
${fitsdir}/lalapps_fits_table_list "WeaveOutNoSim.fits[per_seg_info]" > per_seg_info_no_sim
${fitsdir}/lalapps_fits_table_list "WeaveOutSimFull.fits[per_seg_info]" > per_seg_info_sim_full
diff per_seg_info_no_sim per_seg_info_sim_full
${fitsdir}/lalapps_fits_table_list "WeaveOutNoSim.fits[per_seg_info][col c1=segment_start_s; c2=segment_start_ns; c3=segment_end_s; c4=segment_end_ns; c5=coh_total; c6=coh_total_recomp]" > per_seg_info_no_sim
${fitsdir}/lalapps_fits_table_list "WeaveOutSimMin.fits[per_seg_info][col c1=segment_start_s; c2=segment_start_ns; c3=segment_end_s; c4=segment_end_ns; c5=coh_total; c6=coh_total_recomp]" > per_seg_info_sim_min
diff per_seg_info_no_sim per_seg_info_sim_min
set +x
echo
