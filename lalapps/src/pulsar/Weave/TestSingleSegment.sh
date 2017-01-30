# Perform a fully-coherent search of a single segment, and compare F-statistics to lalapps_ComputeFstatistic_v2

echo "=== Create single-segment search setup ==="
set -x
${builddir}/lalapps_WeaveSetup --first-segment=1122332211/90000 --detectors=H1,L1 --output-file=WeaveSetup.fits
set +x
echo

echo "=== Extract reference time from WeaveSetup.fits ==="
set -x
${fitsdir}/lalapps_fits_header_getval "WeaveSetup.fits[0]" 'DATE-OBS GPS' > tmp
ref_time=`cat tmp | xargs printf "%.9f"`
set +x
echo

echo "=== Generate SFTs with injected signal ==="
set -x
inject_params="Alpha=5.4; Delta=1.1; Freq=70.5; f1dot=-1e-9"
${injdir}/lalapps_Makefakedata_v5 --randSeed=1234 --fmin=70.0 --Band=1.0 \
    --injectionSources="{refTime=${ref_time}; h0=0.5; cosi=0.2; psi=0.4; phi0=0.1; ${inject_params}}" \
    --Tsft=1800 --outSingleSFT --outSFTdir=. --IFOs=H1,L1 --sqrtSX=1,1 \
    --timestampsFiles=${srcdir}/timestamps-1.txt,${srcdir}/timestamps-2.txt
set +x
echo

echo "=== Perform fully-coherent search ==="
set -x
${builddir}/lalapps_Weave --output-file=WeaveOut.fits \
    --toplists=mean2F --toplist-limit=0 --per-detector --misc-info \
    --setup-file=WeaveSetup.fits --sft-files='*.sft' --Fstat-method=DemodBest \
    --freq=70.5~1e-4 --f1dot=-2e-9,0 --semi-max-mismatch=11
set +x
echo

echo "=== Check approximate/computed number of semicoherent templates"
set -x
${fitsdir}/lalapps_fits_header_getval "WeaveOut.fits[0]" 'SEMIAPPX' > tmp
semi_ntmpl=`cat tmp | xargs printf "%d"`
${fitsdir}/lalapps_fits_header_getval "WeaveOut.fits[0]" 'SEMICOMP' > tmp
semi_ncomp=`cat tmp | xargs printf "%d"`
expr ${semi_ncomp} '=' ${semi_ntmpl}
set +x
echo

echo "=== Check average number of semicoherent templates per dimension is more than one"
set -x
for dim in SSKYA SSKYB NU0DOT NU1DOT; do
    ${fitsdir}/lalapps_fits_header_getval "WeaveOut.fits[0]" "SEMIAVG ${dim}" > tmp
    semi_avg_ntmpl_dim=`cat tmp | xargs printf "%d"`
    expr ${semi_avg_ntmpl_dim} '>' 1
done
set +x
echo

echo "=== Check computed number of coherent results"
set -x
${fitsdir}/lalapps_fits_header_getval "WeaveOut.fits[0]" 'TCOHCOMP' > tmp
tot_coh_ncomp=`cat tmp | xargs printf "%d"`
expr ${tot_coh_ncomp} '=' ${semi_ntmpl}
set +x
echo

echo "=== Check that no results were recomputed ==="
set -x
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[per_seg_info][col coh_nrecomp]" > tmp
coh_nrecomp=`cat tmp | sed "/^#/d" | xargs printf "%d"`
expr ${coh_nrecomp} '=' 0
set +x
echo

echo "=== Extract segment start/end times from WeaveSetup.fits ==="
set -x
${fitsdir}/lalapps_fits_table_list "WeaveSetup.fits[segments][col start_s; start_ns][#row == 1]" > tmp
start_time=`cat tmp | sed "/^#/d" | xargs printf "%d.%09d"`
${fitsdir}/lalapps_fits_table_list "WeaveSetup.fits[segments][col end_s; end_ns][#row == 1]" > tmp
end_time=`cat tmp | sed "/^#/d" | xargs printf "%d.%09d"`
set +x
echo

echo "=== Extract template bank and F-statistics from WeaveOut.fits as ASCII table ==="
set -x
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean2F_toplist][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0]" > WeaveBank.txt
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean2F_toplist][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0; c7=DEFNULL(mean2F,-999)]" > WeaveFstats.txt
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean2F_toplist][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0; c7=DEFNULL(mean2F_H1,-999)]" > WeaveFstatsH1.txt
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean2F_toplist][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0; c7=DEFNULL(mean2F_L1,-999)]" > WeaveFstatsL1.txt
set +x
echo

echo "=== Recompute lalapps_Weave F-statistics using lalapps_ComputeFstatistic_v2 ==="
set -x
${fstatdir}/lalapps_ComputeFstatistic_v2 --outputFstat=CFSv2Fstats.txt --outputSingleFstats --refTime=${ref_time} \
    --minStartTime=${start_time} --maxStartTime=${end_time} --DataFiles='*.sft' \
    --TwoFthreshold=0 --FstatMethod=DemodBest --gridType=6 --gridFile=WeaveBank.txt
sed -i '/^%/d' CFSv2Fstats.txt
awk '{print $1, $2, $3, $4, $5, $6, $8}' CFSv2Fstats.txt > CFSv2FstatsH1.txt
awk '{print $1, $2, $3, $4, $5, $6, $9}' CFSv2Fstats.txt > CFSv2FstatsL1.txt
set +x
echo

echo "=== Compare F-statistics from lalapps_Weave to lalapps_ComputeFstatistic_v2 ==="
set -x
${fstatdir}/lalapps_compareFstats --Fname1=WeaveFstats.txt --Fname2=CFSv2Fstats.txt
${fstatdir}/lalapps_compareFstats --Fname1=WeaveFstatsH1.txt --Fname2=CFSv2FstatsH1.txt
${fstatdir}/lalapps_compareFstats --Fname1=WeaveFstatsL1.txt --Fname2=CFSv2FstatsL1.txt
set +x
echo

echo "=== Compute F-statistic at exact injected signal parameters using lalapps_ComputeFstatistic_v2 ==="
set -x
${fstatdir}/lalapps_ComputeFstatistic_v2 --outputFstat=CFSv2Exact.txt --outputSingleFstats --refTime=${ref_time} \
    --minStartTime=${start_time} --maxStartTime=${end_time} --DataFiles='*.sft' \
    --TwoFthreshold=0 --FstatMethod=ResampBest `echo "${inject_params}" | sed 's/^/--/;s/; / --/g'`
set +x
echo

echo "=== Compare F-statistic at exact injected signal parameters with loudest F-statistic found by lalapps_Weave ==="
set -x
coh2F_exact=`cat CFSv2Exact.txt | sed -n '/^[^%]/{p;q}' | awk '{print $7}'`
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean2F_toplist][col c1=mean2F][#row == 1]" > tmp
coh2F_loud=`cat tmp | sed "/^#/d" | xargs printf "%g"`
# Value of 'mean_mu' was calculated by:
#   octapps_run WeaveFstatMismatch --setup-file=TestSingleSegment.testdir/WeaveSetup.fits --spindowns=1 --semi-max-mismatch=11 --coh-max-mismatch=0 --output=mean
mean_mu=0.6102
awk "BEGIN { print mu = ( ${coh2F_exact} - ${coh2F_loud} ) / ${coh2F_exact}; exit ( mu < ${mean_mu} ? 0 : 1 ) }"
set +x
echo
