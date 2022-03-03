# Perform search with mean-2F histogram, and check for consistent results with toplist

export LAL_FSTAT_FFT_PLAN_MODE=ESTIMATE

echo "=== Create search setup with 3 segments ==="
set -x
lalapps_WeaveSetup --first-segment=1122332211/90000 --segment-count=3 --detectors=H1,L1 --output-file=WeaveSetup.fits
lalapps_fits_overview WeaveSetup.fits
set +x
echo

echo "=== Restrict timestamps to segment list in WeaveSetup.fits ==="
set -x
lalapps_fits_table_list 'WeaveSetup.fits[segments][col c1=start_s; col2=end_s]' \
    | awk 'BEGIN { print "/^#/ { print }" } /^#/ { next } { printf "%i <= $1 && $1 <= %i { print }\n", $1, $2 + 1 }' > timestamp-filter.awk
awk -f timestamp-filter.awk all-timestamps-1.txt > timestamps-1.txt
awk -f timestamp-filter.awk all-timestamps-2.txt > timestamps-2.txt
set +x
echo

echo "=== Generate SFTs ==="
set -x
lalapps_Makefakedata_v5 --randSeed=3456 --fmin=49.5 --Band=2.0 --Tsft=1800 \
    --outSingleSFT --outSFTdir=. --IFOs=H1,L1 --sqrtSX=1,1 \
    --timestampsFiles=timestamps-1.txt,timestamps-2.txt
set +x
echo

echo "=== Perform search with mean-2F histogram ==="
set -x
lalapps_Weave --output-file=WeaveOut.fits \
    --toplists=mean2F --toplist-limit=0 --mean2F-hgrm --setup-file=WeaveSetup.fits --sft-files='*.sft' \
    --sky-patch-count=24 --sky-patch-index=0 --freq=50/0.005 --f1dot=-1e-9,0 --semi-max-mismatch=5 --coh-max-mismatch=0.4
lalapps_fits_overview WeaveOut.fits
set +x
echo

echo "=== Compare histogram to toplist ==="
set -x
lalapps_fits_table_list 'WeaveOut.fits[mean2F_hgrm]' | awk '/^#/ { next } { printf "%8.5f %8.5f %4i\n", $1, $2, $3}' > mean2F_hgrm.txt
lalapps_fits_table_list 'WeaveOut.fits[mean2F_toplist]' | awk '/^#/ { next } { ++hgrm[int(10 * $5)] } END { for (j in hgrm) { printf "%8.5f %8.5f %4i\n", 0.1 * j, 0.1 * (j + 1), hgrm[j] } }' | sort -n > mean2F_hgrm_from_toplist.txt
paste mean2F_hgrm.txt mean2F_hgrm_from_toplist.txt | awk '{ if ($1 != $4 || $2 != $5 || $3 - $6 > 3 || $6 - $3 > 3) { exit(1) } }'
set +x
echo
set +e
diff mean2F_hgrm.txt mean2F_hgrm_from_toplist.txt
set -e
