## create temporary directories
mkdir -p broadband/ narrowband/

## run lalapps_Makefakedata_v5 to create fake broadband SFTs
start=1257800000
span=9000
cmdline="lalapps_Makefakedata_v5 --IFOs H1 --sqrtSX 1e-24 --startTime ${start} --duration ${span} --fmin 10 --Band 85 --injectionSources '{Alpha=0.1; Delta=0.4; Freq=30; f1dot=1e-10; h0=1e-24; cosi=0.7; refTime=1257800000}' --outSingleSFT=no --outSFTdir broadband/"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi

## run lalapps_splitSFTs to create narrowband SFTs
for sft in broadband/*.sft; do
    cmdline="lalapps_splitSFTs -fs 10 -fe 95 -fb 8 -n narrowband/ -i $sft"
    if ! eval "$cmdline"; then
        echo "ERROR: something failed when running '$cmdline'"
        exit 1
    fi
done

## check narrowband SFT names
for freq in 10 18 26 34 42 50 58 66 74 82 90; do
    if [ $freq = 90 ]; then
        binwidth="0005Hz1"
    else
        binwidth="0008Hz0"
    fi
    sft="narrowband/H-5_H1_1800SFT_NB_F00${freq}Hz0_W${binwidth}-${start}-${span}.sft"
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
