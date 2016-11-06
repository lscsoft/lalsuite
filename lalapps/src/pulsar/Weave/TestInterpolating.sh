# Perform an interpolating search, and compare F-statistics to lalapps_ComputeFstatistic_v2

echo "=== Generate SFTs ==="
set -x
${injdir}/lalapps_Makefakedata_v5 --randSeed=3456 --fmin=50 --Band=1.0 --Tsft=1800 \
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

echo "=== Perform interpolating search ==="
set -x
${builddir}/lalapps_Weave --output-file=WeaveOut.fits \
    --output-toplist-limit=5000 --output-per-detector --output-per-segment --output-misc-info \
    --setup-file=WeaveSetup.fits --sft-files='*.sft' --Fstat-method=DemodBest \
    --alpha=2.3/0.9 --delta=-1.2/2.3 --freq=50.5/1e-4 --f1dot=-1.5e-9,0 --semi-max-mismatch=0.6 --coh-max-mismatch=0.3
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

echo "=== Extract reference time from WeaveSetup.fits ==="
set -x
${fitsdir}/lalapps_fits_header_getval "WeaveSetup.fits[0]" 'DATE-OBS GPS' > tmp
ref_time=`cat tmp | xargs printf "%.9f"`
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

    echo "=== Segment #${seg}: Extract coherent template bank from WeaveOut.fits as ASCII table ==="
    set -x
    ${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean_twoF_toplist][col c1=seg${seg}_freq; c2=seg${seg}_alpha; c3=seg${seg}_delta; c4=seg${seg}_f1dot; c5=0; c6=0]" > WeaveSeg${seg}CohBank.txt
    set +x
    echo

    echo "=== Segment #${seg}: Extract coherent template bank and coherent F-statistics from WeaveOut.fits as ASCII table ==="
    set -x
    ${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean_twoF_toplist][col c1=seg${seg}_freq; c2=seg${seg}_alpha; c3=seg${seg}_delta; c4=seg${seg}_f1dot; c5=0; c6=0; c7=DEFNULL(seg${seg}_twoF,-999)]" > WeaveSeg${seg}CohBankCohFstats.txt
    ${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean_twoF_toplist][col c1=seg${seg}_freq; c2=seg${seg}_alpha; c3=seg${seg}_delta; c4=seg${seg}_f1dot; c5=0; c6=0; c7=DEFNULL(seg${seg}_twoF_H1,-999)]" > WeaveSeg${seg}CohBankCohFstatsH1.txt
    ${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean_twoF_toplist][col c1=seg${seg}_freq; c2=seg${seg}_alpha; c3=seg${seg}_delta; c4=seg${seg}_f1dot; c5=0; c6=0; c7=DEFNULL(seg${seg}_twoF_L1,-999)]" > WeaveSeg${seg}CohBankCohFstatsL1.txt
    set +x
    echo

    echo "=== Segment #${seg}: Recompute coherent F-statistics using lalapps_ComputeFstatistic_v2 ==="
    set -x
    ${fstatdir}/lalapps_ComputeFstatistic_v2 --outputFstat=CFSv2Seg${seg}Fstats.txt --outputSingleFstats --refTime=${ref_time} --minStartTime=${start_time} --maxStartTime=${end_time} --DataFiles='*.sft' --TwoFthreshold=0 --FstatMethod=DemodBest --gridType=6 --gridFile=WeaveSeg${seg}CohBank.txt
    sed -i '/^%/d' CFSv2Seg${seg}Fstats.txt
    if [ ${seg} -eq 3 ]; then
        awk '{print $1, $2, $3, $4, $5, $6, $7, $8, -999}' CFSv2Seg${seg}Fstats.txt > tmp
        mv -f tmp CFSv2Seg${seg}Fstats.txt
    fi
    awk '{print $1, $2, $3, $4, $5, $6, $8}' CFSv2Seg${seg}Fstats.txt > CFSv2Seg${seg}FstatsH1.txt
    awk '{print $1, $2, $3, $4, $5, $6, $9}' CFSv2Seg${seg}Fstats.txt > CFSv2Seg${seg}FstatsL1.txt
    set +x
    echo

    echo "=== Segment #${seg}: Compare coherent F-statistics from lalapps_Weave to lalapps_ComputeFstatistic_v2 ==="
    set -x
    ${fstatdir}/lalapps_compareFstats --Fname1=WeaveSeg${seg}CohBankCohFstats.txt --Fname2=CFSv2Seg${seg}Fstats.txt
    ${fstatdir}/lalapps_compareFstats --Fname1=WeaveSeg${seg}CohBankCohFstatsH1.txt --Fname2=CFSv2Seg${seg}FstatsH1.txt
    ${fstatdir}/lalapps_compareFstats --Fname1=WeaveSeg${seg}CohBankCohFstatsL1.txt --Fname2=CFSv2Seg${seg}FstatsL1.txt
    set +x
    echo

done

echo "=== Extract semicoherent template bank from WeaveOut.fits as ASCII table ==="
set -x
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean_twoF_toplist][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0]" > WeaveSemiBank.txt
sed -i '/^#/d' WeaveSemiBank.txt
set +x
echo

echo "=== Extract semicoherent F-statistics from WeaveOut.fits as ASCII table ==="
set -x
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean_twoF_toplist][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0; c7=DEFNULL(mean_twoF,-999)]" > WeaveSemiFstats.txt
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean_twoF_toplist][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0; c7=DEFNULL(mean_twoF_H1,-999)]" > WeaveSemiFstatsH1.txt
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[mean_twoF_toplist][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0; c7=DEFNULL(mean_twoF_L1,-999)]" > WeaveSemiFstatsL1.txt
set +x
echo

echo "=== Add coherent F-statistics computed by lalapps_ComputeFstatistic_v2 over segments ==="
set -x
paste WeaveSemiBank.txt CFSv2Seg1Fstats.txt CFSv2Seg2Fstats.txt CFSv2Seg3Fstats.txt > WeaveSemiBankCFSv2AllSegFstats.txt
awk '{print $1, $2, $3, $4, $5, $6, ($13 + $22 + $31) / 3}' WeaveSemiBankCFSv2AllSegFstats.txt > CFSv2SemiFstats.txt
awk '{print $1, $2, $3, $4, $5, $6, ($14 + $23 + $32) / 3}' WeaveSemiBankCFSv2AllSegFstats.txt > CFSv2SemiFstatsH1.txt
awk '{print $1, $2, $3, $4, $5, $6, ($15 + $24      ) / 2}' WeaveSemiBankCFSv2AllSegFstats.txt > CFSv2SemiFstatsL1.txt
set +x
echo

echo "=== Compare semicoherent F-statistics from lalapps_Weave to lalapps_ComputeFstatistic_v2 ==="
set -x
${fstatdir}/lalapps_compareFstats --Fname1=WeaveSemiFstats.txt --Fname2=CFSv2SemiFstats.txt
${fstatdir}/lalapps_compareFstats --Fname1=WeaveSemiFstatsH1.txt --Fname2=CFSv2SemiFstatsH1.txt
${fstatdir}/lalapps_compareFstats --Fname1=WeaveSemiFstatsL1.txt --Fname2=CFSv2SemiFstatsL1.txt
set +x
echo
