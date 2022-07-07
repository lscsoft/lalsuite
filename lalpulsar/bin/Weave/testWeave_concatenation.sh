# Perform several searches, concatenate output files with lalpulsar_WeaveConcat, and check for consistent results

export LAL_FSTAT_FFT_PLAN_MODE=ESTIMATE

echo "=== Create search setup with 2 segments ==="
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

echo "=== Perform searches with mean-2F histogram ==="
set -x
lalpulsar_Weave --output-file=WeaveOut1.fits \
    --toplists=mean2F --toplist-limit=0 --mean2F-hgrm --setup-file=WeaveSetup.fits --sft-files='*.sft' \
    --sky-patch-count=24 --sky-patch-index=0 --freq=50/0.001 --f1dot=-1e-9,0 --semi-max-mismatch=1.2 --coh-max-mismatch=0.4
lalpulsar_fits_overview WeaveOut1.fits
lalpulsar_Weave --output-file=WeaveOut2.fits \
    --toplists=mean2F --toplist-limit=0 --mean2F-hgrm --setup-file=WeaveSetup.fits --sft-files='*.sft' \
    --sky-patch-count=24 --sky-patch-index=8 --freq=50/0.001 --f1dot=-1e-9,0 --semi-max-mismatch=1.2 --coh-max-mismatch=0.4
lalpulsar_fits_overview WeaveOut2.fits
lalpulsar_Weave --output-file=WeaveOut3.fits \
    --toplists=mean2F --toplist-limit=0 --mean2F-hgrm --setup-file=WeaveSetup.fits --sft-files='*.sft' \
    --sky-patch-count=24 --sky-patch-index=22 --freq=50/0.001 --f1dot=-1e-9,0 --semi-max-mismatch=1.2 --coh-max-mismatch=0.4
lalpulsar_fits_overview WeaveOut3.fits
set +x
echo

echo "=== Concatenate output files ==="
set -x
lalpulsar_WeaveConcat --output-result-file=WeaveOutC.fits --input-result-files=WeaveOut1.fits,WeaveOut2.fits,WeaveOut3.fits
lalpulsar_fits_overview WeaveOutC.fits
set +x
echo

echo "=== Check for consistent accumulated fields ==="
set -x
for field in NCOHRES NCOHTPL NSEMITPL; do
    val1=`lalpulsar_fits_header_getval "WeaveOut1.fits" "${field}" | tr '\n\r' '  ' | awk 'NF == 1 {printf "%d", $1}'`
    val2=`lalpulsar_fits_header_getval "WeaveOut2.fits" "${field}" | tr '\n\r' '  ' | awk 'NF == 1 {printf "%d", $1}'`
    val3=`lalpulsar_fits_header_getval "WeaveOut3.fits" "${field}" | tr '\n\r' '  ' | awk 'NF == 1 {printf "%d", $1}'`
    valC=`lalpulsar_fits_header_getval "WeaveOutC.fits" "${field}" | tr '\n\r' '  ' | awk 'NF == 1 {printf "%d", $1}'`
    echo "${val1} ${val2} ${val3} ${valC}" | awk '{ print err = $1 + $2 + $3 - $4; exit ( err == 0 ? 0 : 1 ) }'
done
for field in 'WALL TOTAL' 'CPU TOTAL'; do
    val1=`lalpulsar_fits_header_getval "WeaveOut1.fits" "${field}" | tr '\n\r' '  ' | awk 'NF == 1 {printf "%0.8f", $1}'`
    val2=`lalpulsar_fits_header_getval "WeaveOut2.fits" "${field}" | tr '\n\r' '  ' | awk 'NF == 1 {printf "%0.8f", $1}'`
    val3=`lalpulsar_fits_header_getval "WeaveOut3.fits" "${field}" | tr '\n\r' '  ' | awk 'NF == 1 {printf "%0.8f", $1}'`
    valC=`lalpulsar_fits_header_getval "WeaveOutC.fits" "${field}" | tr '\n\r' '  ' | awk 'NF == 1 {printf "%0.8f", $1}'`
    echo "${val1} ${val2} ${val3} ${valC}" | awk '{ print err = sqrt( ( $1 + $2 + $3 - $4 )^2 ); exit ( err < 1e-5 ? 0 : 1 ) }'
done
set +x
echo


echo "=== Check for consistent concatenated toplists ==="
set -x
lalpulsar_fits_table_list 'WeaveOut1.fits[mean2F_toplist]' | grep -v '^#' > mean2F_toplist_1.txt
lalpulsar_fits_table_list 'WeaveOut2.fits[mean2F_toplist]' | grep -v '^#' > mean2F_toplist_2.txt
lalpulsar_fits_table_list 'WeaveOut3.fits[mean2F_toplist]' | grep -v '^#' > mean2F_toplist_3.txt
lalpulsar_fits_table_list 'WeaveOutC.fits[mean2F_toplist]' | grep -v '^#' > mean2F_toplist_C.txt
cat mean2F_toplist_1.txt mean2F_toplist_2.txt mean2F_toplist_3.txt | sort -n > mean_toplist_123.txt
cat mean2F_toplist_C.txt | sort -n > mean_toplist_CCC.txt
diff mean_toplist_123.txt mean_toplist_CCC.txt
set +x
echo

echo "=== Check for consistent concatenated histograms ==="
set -x
lalpulsar_fits_table_list 'WeaveOut1.fits[mean2F_hgrm]' | grep -v '^#' > mean2F_hgrm_1.txt
lalpulsar_fits_table_list 'WeaveOut2.fits[mean2F_hgrm]' | grep -v '^#' > mean2F_hgrm_2.txt
lalpulsar_fits_table_list 'WeaveOut3.fits[mean2F_hgrm]' | grep -v '^#' > mean2F_hgrm_3.txt
lalpulsar_fits_table_list 'WeaveOutC.fits[mean2F_hgrm]' | grep -v '^#' > mean2F_hgrm_C.txt
cat mean2F_hgrm_1.txt mean2F_hgrm_2.txt mean2F_hgrm_3.txt | awk '{ hgrm[$1] += $3 } END { for (x in hgrm) { printf "%8.5f %8.5f %4i\n", x, x + 0.1, hgrm[x] } }' | sort -n > mean_hgrm_123.txt
cat mean2F_hgrm_C.txt | awk '{ hgrm[$1] += $3 } END { for (x in hgrm) { printf "%8.5f %8.5f %4i\n", x, x + 0.1, hgrm[x] } }' | sort -n > mean_hgrm_CCC.txt
diff mean_hgrm_123.txt mean_hgrm_CCC.txt
set +x
echo
