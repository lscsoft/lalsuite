# Perform a non-interpolating search, and compare F-statistics to lalapps_ComputeFstatistic_v2

echo "=== Create search setup with 3 segments ==="
set -x
${builddir}/lalapps_WeaveSetup --ref-time=1122334444 --first-segment=1122332211/90000 --segment-count=3 --segment-gap=30000 --detectors=H1,L1 --output-file=WeaveSetup.fits
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
inject_params="Alpha=2.9; Delta=0.71; Freq=50.5; f1dot=-3e-10"
${injdir}/lalapps_Makefakedata_v5 --randSeed=2345 --fmin=50.0 --Band=1.0 \
    --injectionSources="{refTime=${ref_time}; h0=0.5; cosi=0.1; psi=4.4; phi0=2.1; ${inject_params}}" \
    --Tsft=1800 --outSingleSFT --outSFTdir=. --IFOs=H1,L1 --sqrtSX=1,1 \
    --timestampsFiles=${srcdir}/timestamps-regular.txt,${srcdir}/timestamps-irregular.txt
set +x
echo

echo "=== Perform non-interpolating search ==="
set -x
${builddir}/lalapps_Weave --output-file=WeaveOut.fits \
    --toplist-limit=3000 --per-detector --per-segment --misc-info \
    --setup-file=WeaveSetup.fits --sft-files='*.sft' --Fstat-method=DemodBest \
    --alpha=2.3/0.9 --delta=-1.2/2.3 --freq=50.5~1e-4 --f1dot=-1e-9,0 --semi-max-mismatch=6 --interpolation=no
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
${fitsdir}/lalapps_fits_header_getval "WeaveOut.fits[0]" 'NSEGMENT' > tmp
nsegments=`cat tmp | xargs printf "%d"`
${fitsdir}/lalapps_fits_header_getval "WeaveOut.fits[0]" 'TCOHCOMP' > tmp
tot_coh_ncomp=`cat tmp | xargs printf "%d"`
expr ${tot_coh_ncomp} '=' ${semi_ntmpl} '*' ${nsegments}
set +x
echo

echo "=== Extract semicoherent template bank from WeaveOut.fits as ASCII table ==="
set -x
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean2F_toplist][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0]" > WeaveSemiBank.txt
set +x
echo

for seg in 1 2 3; do

    echo "=== Segment #${seg}: Check that no results were recomputed ==="
    set -x
    ${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[per_seg_info][col coh_nrecomp][#row == ${seg}]" > tmp
    coh_nrecomp=`cat tmp | sed "/^#/d" | xargs printf "%d"`
    expr ${coh_nrecomp} '=' 0
    set +x
    echo

    echo "=== Segment #${seg}: Extract segment start/end times from WeaveSetup.fits ==="
    set -x
    ${fitsdir}/lalapps_fits_table_list "WeaveSetup.fits[segments][col start_s; start_ns][#row == ${seg}]" > tmp
    start_time=`cat tmp | sed "/^#/d" | xargs printf "%d.%09d"`
    ${fitsdir}/lalapps_fits_table_list "WeaveSetup.fits[segments][col end_s; end_ns][#row == ${seg}]" > tmp
    end_time=`cat tmp | sed "/^#/d" | xargs printf "%d.%09d"`
    set +x
    echo

    echo "=== Segment #${seg}: Extract semicoherent template bank and coherent F-statistics from WeaveOut.fits as ASCII table ==="
    set -x
    ${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean2F_toplist][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0; c7=DEFNULL(coh2F_seg[${seg}],-999)]" > WeaveSeg${seg}SemiBankCohFstats.txt
    ${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean2F_toplist][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0; c7=DEFNULL(coh2F_H1_seg[${seg}],-999)]" > WeaveSeg${seg}SemiBankCohFstatsH1.txt
    ${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean2F_toplist][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0; c7=DEFNULL(coh2F_L1_seg[${seg}],-999)]" > WeaveSeg${seg}SemiBankCohFstatsL1.txt
    set +x
    echo

    echo "=== Segment #${seg}: Extract coherent template bank and coherent F-statistics from WeaveOut.fits as ASCII table ==="
    set -x
    ${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean2F_toplist][col c1=freq_seg[${seg}]; c2=alpha_seg[${seg}]; c3=delta_seg[${seg}]; c4=f1dot_seg[${seg}]; c5=0; c6=0; c7=DEFNULL(coh2F_seg[${seg}],-999)]" > WeaveSeg${seg}CohBankCohFstats.txt
    ${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean2F_toplist][col c1=freq_seg[${seg}]; c2=alpha_seg[${seg}]; c3=delta_seg[${seg}]; c4=f1dot_seg[${seg}]; c5=0; c6=0; c7=DEFNULL(coh2F_H1_seg[${seg}],-999)]" > WeaveSeg${seg}CohBankCohFstatsH1.txt
    ${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean2F_toplist][col c1=freq_seg[${seg}]; c2=alpha_seg[${seg}]; c3=delta_seg[${seg}]; c4=f1dot_seg[${seg}]; c5=0; c6=0; c7=DEFNULL(coh2F_L1_seg[${seg}],-999)]" > WeaveSeg${seg}CohBankCohFstatsL1.txt
    set +x
    echo

    echo "=== Segment #${seg}: Check that semicoherent and coherent template banks agree ==="
    set -x
    ${fstatdir}/lalapps_compareFstats --Fname1=WeaveSeg${seg}SemiBankCohFstats.txt --Fname2=WeaveSeg${seg}CohBankCohFstats.txt
    ${fstatdir}/lalapps_compareFstats --Fname1=WeaveSeg${seg}SemiBankCohFstatsH1.txt --Fname2=WeaveSeg${seg}CohBankCohFstatsH1.txt
    ${fstatdir}/lalapps_compareFstats --Fname1=WeaveSeg${seg}SemiBankCohFstatsL1.txt --Fname2=WeaveSeg${seg}CohBankCohFstatsL1.txt
    set +x
    echo

    echo "=== Segment #${seg}: Recompute coherent F-statistics using lalapps_ComputeFstatistic_v2 ==="
    set -x
    ${fstatdir}/lalapps_ComputeFstatistic_v2 --outputFstat=CFSv2Seg${seg}Fstats.txt --outputSingleFstats --refTime=${ref_time} \
        --minStartTime=${start_time} --maxStartTime=${end_time} --DataFiles='*.sft' \
        --TwoFthreshold=0 --FstatMethod=DemodBest --gridType=6 --gridFile=WeaveSemiBank.txt
    sed -i '/^%/d' CFSv2Seg${seg}Fstats.txt
    if [ ${seg} -eq 3 ]; then
        awk '{print $1, $2, $3, $4, $5, $6, $7, -999, $8}' CFSv2Seg${seg}Fstats.txt > tmp
        mv -f tmp CFSv2Seg${seg}Fstats.txt
    fi
    awk '{print $1, $2, $3, $4, $5, $6, $8}' CFSv2Seg${seg}Fstats.txt > CFSv2Seg${seg}FstatsH1.txt
    awk '{print $1, $2, $3, $4, $5, $6, $9}' CFSv2Seg${seg}Fstats.txt > CFSv2Seg${seg}FstatsL1.txt
    set +x
    echo

    echo "=== Segment #${seg}: Compare coherent F-statistics from lalapps_Weave to lalapps_ComputeFstatistic_v2 ==="
    set -x
    ${fstatdir}/lalapps_compareFstats --Fname1=WeaveSeg${seg}SemiBankCohFstats.txt --Fname2=CFSv2Seg${seg}Fstats.txt
    ${fstatdir}/lalapps_compareFstats --Fname1=WeaveSeg${seg}SemiBankCohFstatsH1.txt --Fname2=CFSv2Seg${seg}FstatsH1.txt
    ${fstatdir}/lalapps_compareFstats --Fname1=WeaveSeg${seg}SemiBankCohFstatsL1.txt --Fname2=CFSv2Seg${seg}FstatsL1.txt
    set +x
    echo

    echo "=== Segment #${seg}: Compute F-statistic at exact injected signal parameters using lalapps_ComputeFstatistic_v2 ==="
    set -x
    ${fstatdir}/lalapps_ComputeFstatistic_v2 --outputFstat=CFSv2Seg${seg}Exact.txt --outputSingleFstats --refTime=${ref_time} \
        --minStartTime=${start_time} --maxStartTime=${end_time} --DataFiles='*.sft' \
        --TwoFthreshold=0 --FstatMethod=ResampBest `echo "${inject_params}" | sed 's/^/--/;s/; / --/g'`
    set +x
    echo

done

echo "=== Extract semicoherent F-statistics from WeaveOut.fits as ASCII table ==="
set -x
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean2F_toplist][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0; c7=DEFNULL(mean2F,-999)]" > WeaveSemiFstats.txt
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean2F_toplist][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0; c7=DEFNULL(mean2F_H1,-999)]" > WeaveSemiFstatsH1.txt
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean2F_toplist][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0; c7=DEFNULL(mean2F_L1,-999)]" > WeaveSemiFstatsL1.txt
set +x
echo

echo "=== Add coherent F-statistics computed by lalapps_ComputeFstatistic_v2 over segments ==="
set -x
paste CFSv2Seg1Fstats.txt CFSv2Seg2Fstats.txt CFSv2Seg3Fstats.txt > CFSv2AllSegFstats.txt
awk '{print $1, $2, $3, $4, $5, $6, ($7 + $16 + $25) / 3}' CFSv2AllSegFstats.txt > CFSv2SemiFstats.txt
awk '{print $1, $2, $3, $4, $5, $6, ($8 + $17      ) / 2}' CFSv2AllSegFstats.txt > CFSv2SemiFstatsH1.txt
awk '{print $1, $2, $3, $4, $5, $6, ($9 + $18 + $27) / 3}' CFSv2AllSegFstats.txt > CFSv2SemiFstatsL1.txt
set +x
echo

echo "=== Compare semicoherent F-statistics from lalapps_Weave to lalapps_ComputeFstatistic_v2 ==="
set -x
${fstatdir}/lalapps_compareFstats --Fname1=WeaveSemiFstats.txt --Fname2=CFSv2SemiFstats.txt
${fstatdir}/lalapps_compareFstats --Fname1=WeaveSemiFstatsH1.txt --Fname2=CFSv2SemiFstatsH1.txt
${fstatdir}/lalapps_compareFstats --Fname1=WeaveSemiFstatsL1.txt --Fname2=CFSv2SemiFstatsL1.txt
set +x
echo

echo "=== Compare F-statistic at exact injected signal parameters with loudest F-statistic found by lalapps_Weave ==="
set -x
paste CFSv2Seg1Exact.txt CFSv2Seg2Exact.txt CFSv2Seg3Exact.txt > CFSv2AllSegExact.txt
coh2F_exact=`cat CFSv2AllSegExact.txt | sed -n '/^[^%]/{p;q}' | awk '{print ($7 + $16 + $25) / 3}'`
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean2F_toplist][col c1=mean2F][#row == 1]" > tmp
coh2F_loud=`cat tmp | sed "/^#/d" | xargs printf "%g"`
# Value of 'mean_mu' was calculated by:
#   octapps_run WeaveFstatMismatch --setup-file=TestNonInterpolating.testdir/WeaveSetup.fits --spindowns=1 --semi-max-mismatch=6 --coh-max-mismatch=0 --output=mean
mean_mu=0.52
awk "BEGIN { print mu = ( ${coh2F_exact} - ${coh2F_loud} ) / ${coh2F_exact}; exit ( mu < ${mean_mu} ? 0 : 1 ) }"
set +x
echo
