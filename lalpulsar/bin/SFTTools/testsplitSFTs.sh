## run lalpulsar_Makefakedata_v5 to create fake broadband SFTs
Tsft=1800
start=1257800000
span=9000
misc='SFT_desc_w_us'   # test that 'misc' field can contain underscores
cmdlinebase="lalpulsar_Makefakedata_v5 --outLabel '${misc}' --IFOs H1 --sqrtSX 1e-24 --startTime ${start} --duration ${span} --fmin 10 --Band 85 --injectionSources '{Alpha=0.1; Delta=0.4; Freq=30; f1dot=1e-10; h0=1e-24; cosi=0.7; refTime=1257800000}'"
mkdir -p broadband1a/
cmdline="${cmdlinebase} --outSingleSFT=no --outSFTdir broadband1a/"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
ln -s broadband1a broadband1b
mkdir -p broadband1c/
cmdline="${cmdlinebase} --outSingleSFT=yes --outSFTdir broadband1c/"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi

## run lalpulsar_splitSFTs to create narrowband SFTs
mkdir -p narrowband1a/
cmdline="lalpulsar_splitSFTs -fs 10 -fe 95 -fb 8 -n narrowband1a/"
for sft in broadband1a/*.sft; do
    cmdline="${cmdline} $sft"
done
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
mkdir -p narrowband1b/
for sft in broadband1b/*.sft; do
    cmdline="lalpulsar_splitSFTs -fs 10 -fe 95 -fb 8 -n narrowband1b/ $sft"
    if ! eval "$cmdline"; then
        echo "ERROR: something failed when running '$cmdline'"
        exit 1
    fi
done
mkdir -p narrowband1c/
cmdline="lalpulsar_splitSFTs -fs 10 -fe 95 -fb 8 -n narrowband1c/ -- broadband1c/H-5_H1_${Tsft}SFT_${misc}-${start}-${span}.sft"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi

for i in 1a 1b 1c; do

    ## check narrowband SFT names
    for freq in 10 18 26 34 42 50 58 66 74 82 90; do
        if [ $freq = 90 ]; then
            binwidth="0005Hz0"
        else
            binwidth="0008Hz0"
        fi
        sft="narrowband${i}/H-5_H1_${Tsft}SFT_NBF00${freq}Hz0W${binwidth}_${misc}-${start}-${span}.sft"
        if ! test -f $sft; then
            echo "ERROR: could not find file '$sft'"
            exit 1
        fi
    done

    ## dump broadband and narrowband SFT bins, sorted by frequency
    for sftband in broadband narrowband; do
        cmdline="lalpulsar_dumpSFT -d -i '${sftband}${i}/*' | sed '/^$/d;/^%/d' | sort -n -k1,3 > ${sftband}${i}_bins.txt"
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

## split broadband1c/ set again, but with timestamps constraints

endTime1=$(echo "$start + $span" | bc)
startTime2=$(echo "$start + $Tsft" | bc)
endTime2=$(echo "$endTime1 - $Tsft" | bc)
span2=$(echo "$endTime2 - $startTime2" | bc)

## first using the constraint user options, but trying to get back the full split set
mkdir -p narrowband2a/
cmdline="lalpulsar_splitSFTs -fs 10 -fe 95 -fb 8 -n narrowband2a/ -ts ${start} -te ${endTime1} -- broadband1c/H-5_H1_${Tsft}SFT_${misc}-${start}-${span}.sft"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
## then actually reduce the range
mkdir -p narrowband2b/
cmdline="lalpulsar_splitSFTs -fs 10 -fe 95 -fb 8 -n narrowband2b/ -ts ${startTime2} -te ${endTime2} -- broadband1c/H-5_H1_${Tsft}SFT_${misc}-${start}-${span}.sft"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi

## check narrowband2a/ SFT names (full range)
for freq in 10 18 26 34 42 50 58 66 74 82 90; do
    if [ $freq = 90 ]; then
        binwidth="0005Hz0"
    else
        binwidth="0008Hz0"
    fi
    sft="narrowband2a/H-5_H1_${Tsft}SFT_NBF00${freq}Hz0W${binwidth}_${misc}-${start}-${span}.sft"
    if ! test -f $sft; then
        echo "ERROR: could not find file '$sft'"
        exit 1
    fi
done

## dump narrowband2a/ (full range) SFT bins, sorted by frequency, should match broadband1c/ from before
cmdline="lalpulsar_dumpSFT -d -i 'narrowband2a/*' | sed '/^$/d;/^%/d' | sort -n -k1,3 > narrowband2a_bins.txt"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi

if ! diff -s narrowband1c_bins.txt narrowband2a_bins.txt; then
    echo "ERROR: narrowband1c_bins.txt and narrowband2a_bins.txt should be equal"
    exit 1
fi

## check narrowband2b/ SFT names (reduced range)
for freq in 10 18 26 34 42 50 58 66 74 82 90; do
    if [ $freq = 90 ]; then
        binwidth="0005Hz0"
    else
        binwidth="0008Hz0"
    fi
    sft="narrowband2b/H-3_H1_${Tsft}SFT_NBF00${freq}Hz0W${binwidth}_${misc}-${startTime2}-${span2}.sft"
    if ! test -f $sft; then
        echo "ERROR: could not find file '$sft'"
        exit 1
    fi
done

## check rounding of bins
mkdir -p broadband3/
cmdline="lalpulsar_Makefakedata_v5 --outLabel '${misc}' --IFOs H1 --sqrtSX 1e-24 --startTime ${start} --duration ${Tsft} --Tsft ${Tsft} --fmin 50 --Band 10 --outSingleSFT=no --outSFTdir broadband3/"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
mkdir -p narrowband3a/
cmdline="lalpulsar_splitSFTs --output-directory narrowband3a/ --start-frequency 51 --frequency-bandwidth 1 --end-frequency 52 -- broadband3/H-1_H1_1800SFT_${misc}-1257800000-1800.sft"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
files=`echo narrowband3a/*.sft`
if [ "X$files" != "Xnarrowband3a/H-1_H1_1800SFT_NBF0051Hz0W0001Hz0_${misc}-1257800000-1800.sft" ]; then
    echo "ERROR: extra SFTs generated in narrowband3a/: ${files}"
    exit 1
fi
mkdir -p narrowband3b/
cmdline="lalpulsar_splitSFTs --output-directory narrowband3b/ --start-frequency 51 --frequency-bandwidth 1.2212 --end-frequency 52.2212 -- broadband3/H-1_H1_1800SFT_${misc}-1257800000-1800.sft"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
files=`echo narrowband3b/*.sft`
if [ "X$files" != "Xnarrowband3b/H-1_H1_1800SFT_NBF0051Hz0W0001Hz399_${misc}-1257800000-1800.sft" ]; then
    echo "ERROR: extra SFTs generated in narrowband3b/: ${files}"
    exit 1
fi
mkdir -p narrowband3c/
cmdline="lalpulsar_splitSFTs --output-directory narrowband3c/ --start-frequency 50.9508 --frequency-bandwidth 1.2212 --end-frequency 52.172 -- broadband3/H-1_H1_1800SFT_${misc}-1257800000-1800.sft"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
files=`echo narrowband3c/*.sft`
if [ "X$files" != "Xnarrowband3c/H-1_H1_1800SFT_NBF0050Hz1711W0001Hz399_${misc}-1257800000-1800.sft" ]; then
    echo "ERROR: extra SFTs generated in narrowband3c/: ${files}"
    exit 1
fi
