# Perform a fully-coherent search of a single segment, and compare F-statistics to lalapps_ComputeFstatistic_v2

echo "=== Generate SFTs ==="
set -x
${injdir}/lalapps_Makefakedata_v5 --randSeed=1234 --fmin=100 --Band=1.0 \
    --injectionSources="{refTime=1122332211; h0=0.5; cosi=0.2; psi=0.4; phi0=0.1; Alpha=5.4; Delta=1.1; Freq=100.5; f1dot=-1e-9}" \
    --Tsft=1800 --outSingleSFT --outSFTdir=. --IFOs=H1,L1 --sqrtSX=1,1 \
    --timestampsFiles=${srcdir}/timestamps-irregular.txt,${srcdir}/timestamps-regular.txt
set +x
echo

echo "=== Create single-segment search setup ==="
set -x
${builddir}/lalapps_WeaveSetup --first-segment=1122332211/90000 --detectors=H1,L1 --output-file=WeaveSetup.fits
set +x
echo

echo "=== Perform fully-coherent search ==="
set -x
${builddir}/lalapps_Weave --output-file=WeaveOut.fits --output-max-size=0 --output-per-detector --setup-file=WeaveSetup.fits \
    --sft-files='*.sft' --freq=100.5 --f1dot=-1e-8,0 --semi-max-mismatch=0.3 --Fstat-method=DemodBest
set +x
echo

echo "=== Extract reference time and segment start/end times from WeaveSetup.fits ==="
set -x
${fitsdir}/lalapps_fits_header_list "WeaveSetup.fits[0]" > tmp
ref_time=`cat tmp | sed -n "s/COMMENT DATE-OBS = GPS //p"`
${fitsdir}/lalapps_fits_table_list "WeaveSetup.fits[segments][col start_s; start_ns][#row == 1]" > tmp
start_time=`cat tmp | sed "/^#/d" | xargs printf "%d.%09d"`
${fitsdir}/lalapps_fits_table_list "WeaveSetup.fits[segments][col end_s; end_ns][#row == 1]" > tmp
end_time=`cat tmp | sed "/^#/d" | xargs printf "%d.%09d"`
set +x
echo

echo "=== Extract template bank and F-statistics from WeaveOut.fits as ASCII table ==="
set -x
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[toplist_mean_twoF][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0]" > WeaveBank.txt
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[toplist_mean_twoF][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0; c7=DEFNULL(mean_twoF,-999)]" > WeaveFstats.txt
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[toplist_mean_twoF][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0; c7=DEFNULL(mean_twoF_H1,-999)]" > WeaveFstatsH1.txt
${fitsdir}/lalapps_fits_table_list "WeaveOut.fits[toplist_mean_twoF][col c1=freq; c2=alpha; c3=delta; c4=f1dot; c5=0; c6=0; c7=DEFNULL(mean_twoF_L1,-999)]" > WeaveFstatsL1.txt
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
