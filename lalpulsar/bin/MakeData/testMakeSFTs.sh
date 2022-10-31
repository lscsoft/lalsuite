if test "${LALFRAME_ENABLED}" = false; then
    echo "Skipping test: requires LALFrame"
    exit 77
fi

## common variables
Tsft=1800
Band=256

## make several frames/SFTs with a gap in between
tstart1=1257741529
duration1=$(echo "${Tsft} * 3" | bc)
tend1=$(echo "${tstart1} + ${duration1}" | bc)
tstart2=$(echo "${tend1} + 3215" | bc )
duration2=$(echo "${Tsft} * 2" | bc)
tend2=$(echo "${tstart2} + ${duration2}" | bc)
timestamps="${tstart1}"
for i in 1 2; do
    timestamps="${timestamps} $(echo "${tstart1} + ${Tsft}*${i}" | bc)"
done
for i in 0 1; do
    timestamps="${timestamps} $(echo "${tstart2} + ${Tsft}*${i}" | bc)"
done

## run MFDv5 to create a fake frames and equivalent SFTs
mkdir -p MFDv5/
MFDv5_cmdline_base="lalpulsar_Makefakedata_v5 --IFOs H1 --sqrtSX 1e-24 --Tsft ${Tsft} --fmin 0 --Band ${Band} --injectionSources '{Alpha=0.1; Delta=0.4; Freq=50; f1dot=1e-10; h0=1e-24; cosi=0.7; refTime=${tstart1}}' --outSingleSFT false"
cmdline="${MFDv5_cmdline_base} --outFrameDir MFDv5/ --outSFTdir MFDv5/ --startTime ${tstart1} --duration ${duration1}"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
cmdline="${MFDv5_cmdline_base} --outFrameDir MFDv5/ --outSFTdir MFDv5/ --startTime ${tstart2} --duration ${duration2}"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
MFDv5gwfs="./MFDv5/H-H1_mfdv5-${tstart1}-${duration1}.gwf ./MFDv5/H-H1_mfdv5-${tstart2}-${duration2}.gwf"
for MFDv5gwf in $MFDv5gwfs; do
    if ! test -f $MFDv5gwf; then
        echo "ERROR: could not find file '$MFDv5gwf'"
        exit 1
    fi
done
for ts in ${timestamps}; do
    MFDv5sft="./MFDv5/H-1_H1_${Tsft}SFT_mfdv5-${ts}-${Tsft}.sft"
    if ! test -f $MFDv5sft; then
        echo "ERROR: could not find file '$MFDv5sft'"
        exit 1
    fi
done

## make a cache file for the MFDv5 frame
framecache="./framecache"
echo "H H1_mfdv5 ${tstart1} ${duration1} file://localhost$PWD/MFDv5/H-H1_mfdv5-${tstart1}-${duration1}.gwf"  > $framecache
echo "H H1_mfdv5 ${tstart2} ${duration2} file://localhost$PWD/MFDv5/H-H1_mfdv5-${tstart2}-${duration2}.gwf" >> $framecache

## run MakeSFTs to create SFTs from the fake frames
mkdir -p MSFTs/
pubObsRun=4
pubObsKind="DEV"
pubVersion=1
MSFTs_cmdline_base="lalpulsar_MakeSFTs --frame-cache $framecache --channel-name H1:mfdv5 --sft-duration ${Tsft} --high-pass-freq 0 --start-freq 0 --band ${Band} --comment-field 'Test comment'"
MSFTs_cmdline_public="${MSFTs_cmdline_base} --observing-run ${pubObsRun} --observing-kind ${pubObsKind} --observing-version ${pubVersion}"
cmdline="${MSFTs_cmdline_public} --sft-write-path MSFTs/ --gps-start-time ${tstart1} --gps-end-time ${tend2} --window-type rectangular"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
pubField="O${pubObsRun}${pubObsKind}+V${pubVersion}+Cmfdv5+WRECT"
for ts in ${timestamps}; do
    MSFTsft="./MSFTs/H-1_H1_${Tsft}SFT_${pubField}-${ts}-${Tsft}.sft"
    if ! test -f $MSFTsft; then
        echo "ERROR: could not find file '$MSFTsft'"
        exit 1
    fi
done

## compare SFTs produced by MFDv5 and MakeSFTs, should be very nearly identical
tol=1e-10
for ts in ${timestamps}; do
    MFDv5sft="./MFDv5/H-1_H1_${Tsft}SFT_mfdv5-${ts}-${Tsft}.sft"
    MSFTsft="./MSFTs/H-1_H1_${Tsft}SFT_${pubField}-${ts}-${Tsft}.sft"
    cmdline="lalpulsar_compareSFTs -V -e $tol -1 $MFDv5sft -2 $MSFTsft"
    echo "Comparing SFTs produced by MFDv5 and MakeSFTs, allowed tolerance=$tol:"
    if ! eval "$cmdline"; then
        echo "ERROR: something failed when running '$cmdline'"
        exit 1
    fi
done

## run MakeSFTs to create overlapped SFTs from the 1st fake frame
mkdir -p MSFTs-overlapped/
MSFTs_cmdline_private="${MSFTs_cmdline_base} --observing-run 0 --misc-desc private"
cmdline="${MSFTs_cmdline_private} --sft-write-path MSFTs-overlapped/ --gps-start-time ${tstart1} --gps-end-time ${tend1} --window-type rectangular --overlap-fraction 0.5"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
for i in 0.0 0.5 1.0 1.5 2.0; do
    ts=$(echo "( ${tstart1} + ${Tsft}*${i} ) / 1" | bc)
    MSFTsft="./MSFTs-overlapped/H-1_H1_${Tsft}SFT_private-${ts}-${Tsft}.sft"
    if ! test -f $MSFTsft; then
        echo "ERROR: could not find file '$MSFTsft'"
        exit 1
    fi
done
