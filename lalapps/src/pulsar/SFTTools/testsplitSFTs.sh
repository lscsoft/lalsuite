## create temporary directories
mkdir -p broadband/ narrowband/

## run lalapps_Makefakedata_v5 to create fake broadband SFTs
cmdline="lalapps_Makefakedata_v5 --IFOs H1 --sqrtSX 1e-24 --startTime 1257800000 --duration 9000 --fmin 10 --Band 85 --injectionSources '{Alpha=0.1; Delta=0.4; Freq=30; f1dot=1e-10; h0=1e-24; cosi=0.7; refTime=1257800000}' --outSingleSFT=no --outSFTdir broadband/"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi

## run lalapps_splitSFTs to create narrowband SFTs
for sft in broadband/*.sft; do
    cmdline="lalapps_splitSFTs -fs 10 -fe 95 -fb 8 -o narrowband/ -i $sft"
    if ! eval "$cmdline"; then
        echo "ERROR: something failed when running '$cmdline'"
        exit 1
    fi
done

## check narrowband SFT names
for bin in 18000 32400 46800 61200 75600 90000 104400 118800 133200 147600 162000; do
    sft="narrowband/${bin}"
    if ! test -f $sft; then
        echo "ERROR: could not find file '$sft'"
        exit 1
    fi
done

## dump broadband and narrowband SFT bins, sorted by frequency
for sftband in broadband narrowband; do
    cmdline="lalapps_dumpSFT -d -i '${sftband}/*' | sed '/^$/d;/^%/d' | sort -n -k1,3 > ${sftband}_bins.txt"
    if ! eval "$cmdline"; then
        echo "ERROR: something failed when running '$cmdline'"
        exit 1
    fi
done

## broadband and narrowband SFT bins should be equal
if ! diff -s broadband_bins.txt narrowband_bins.txt; then
    echo "ERROR: broadband_bins.txt and narrowband_bins.txt should be equal"
    exit 1
fi
