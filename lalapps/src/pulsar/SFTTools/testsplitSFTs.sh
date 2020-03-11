## create temporary directories
mkdir -p broadband1/ broadband2/ narrowband1/ narrowband2/

## run lalapps_Makefakedata_v5 to create fake broadband SFTs
start=1257800000
span=9000
cmdlinebase="lalapps_Makefakedata_v5 --IFOs H1 --sqrtSX 1e-24 --startTime ${start} --duration ${span} --fmin 10 --Band 85 --injectionSources '{Alpha=0.1; Delta=0.4; Freq=30; f1dot=1e-10; h0=1e-24; cosi=0.7; refTime=1257800000}'"
cmdline="${cmdlinebase} --outSingleSFT=no --outSFTdir broadband1/"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
cmdline="${cmdlinebase} --outSingleSFT=yes --outSFTdir broadband2/"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi

## run lalapps_splitSFTs to create narrowband SFTs
cmdline="lalapps_splitSFTs -fs 10 -fe 95 -fb 8 -n narrowband1/"
for sft in broadband1/*.sft; do
    cmdline="${cmdline} $sft"
done
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
cmdline="lalapps_splitSFTs -fs 10 -fe 95 -fb 8 -n narrowband2/ -- broadband2/H-5_H1_1800SFT_mfdv5-${start}-${span}.sft"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi

for i in 1 2; do

## check narrowband SFT names
for freq in 10 18 26 34 42 50 58 66 74 82 90; do
    if [ $freq = 90 ]; then
        binwidth="0005Hz1"
    else
        binwidth="0008Hz0"
    fi
    sft="narrowband${i}/H-5_H1_1800SFT_NB_F00${freq}Hz0_W${binwidth}-${start}-${span}.sft"
    if ! test -f $sft; then
        echo "ERROR: could not find file '$sft'"
        exit 1
    fi
done

## dump broadband and narrowband SFT bins, sorted by frequency
for sftband in broadband narrowband; do
    cmdline="lalapps_dumpSFT -d -i '${sftband}${i}/*' | sed '/^$/d;/^%/d' | sort -n -k1,3 > ${sftband}${i}_bins.txt"
    if ! eval "$cmdline"; then
        echo "ERROR: something failed when running '$cmdline'"
        exit 1
    fi
done

## broadband and narrowband SFT bins should be equal
if ! diff -s broadband${i}_bins.txt narrowband${i}_bins.txt; then
    echo "ERROR: broadband${i}_bins.txt and narrowband${i}_bins.txt should be equal"
    exit 1
fi

done
