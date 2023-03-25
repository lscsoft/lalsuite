if test "${LALFRAME_ENABLED}" = false; then
    echo "Skipping test: requires LALFrame"
    exit 77
fi

## common variables
tstart=1257741529
Tsft=1800
tend=$(echo "${tstart} + ${Tsft}" | bc)
Band=256

## run MFDv5 to create fake frame
mkdir -p MFDv5/
cmdline="lalpulsar_Makefakedata_v5 --IFOs H1 --sqrtSX 1e-24 --Tsft ${Tsft} --fmin 0 --Band ${Band} --injectionSources '{Alpha=0.1; Delta=0.4; Freq=50; f1dot=1e-10; h0=1e-24; cosi=0.7; refTime=${tstart}}' --outSingleSFT false --outFrameDir MFDv5/ --startTime ${tstart} --duration ${Tsft}"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
MFDv5gwf="./MFDv5/H-H1_mfdv5-${tstart}-${Tsft}.gwf"
if ! test -f $MFDv5gwf; then
    echo "ERROR: could not find file '$MFDv5gwf'"
    exit 1
fi

## make a cache file for the MFDv5 frame
framecache="./framecache"
echo "H H1_mfdv5 ${tstart} ${Tsft} file://localhost$PWD/MFDv5/H-H1_mfdv5-${tstart}-${Tsft}.gwf" > $framecache

## create multi-channel frames
mkdir -p MultiChFrame/
makeMultiChFrame ./framecache H1:mfdv5 MultiChFrame/
MCTgwf="./MultiChFrame/H-H1_MultiChFrame-${tstart}-${Tsft}.gwf"
if ! test -f $MCFgwf; then
    echo "ERROR: could not find file '$MCFgwf'"
    exit 1
fi

# loop over multiple frame channels
for MCFch in AdcREAL4 AdcREAL8 ProcREAL4 ProcREAL8 SimREAL4 SimREAL8; do
    echo "=== frame channel: ${MCFch} ==="

    ## make a cache file for the MultiChFrame frames
    MCFframecache="./MCFframecache"
    echo "H H1_MultiChFrame ${tstart} 1800 file://localhost$PWD/MultiChFrame/H-H1_MultiChFrame-${tstart}-${Tsft}.gwf" > $MCFframecache

    ## run MakeSFTs to create SFTs from the fake frames
    rm -rf MSFTs/
    mkdir -p MSFTs/
    pubObsRun=4
    pubObsKind="DEV"
    pubRevision=1
    cmdline="lalpulsar_MakeSFTs --frame-cache $MCFframecache --channel-name H1:${MCFch} --sft-duration ${Tsft} --high-pass-freq 0 --start-freq 0 --band ${Band} --comment-field 'Test comment' --observing-run ${pubObsRun} --observing-kind ${pubObsKind} --observing-revision ${pubRevision} --sft-write-path MSFTs/ --gps-start-time ${tstart} --gps-end-time ${tend} --window-type rectangular"
    if ! eval "$cmdline"; then
        echo "ERROR: something failed when running '$cmdline'"
        exit 1
    fi
    pubField="O${pubObsRun}${pubObsKind}+R${pubRevision}+C${MCFch}+WRECT"
    MSFTsft="./MSFTs/H-1_H1_${Tsft}SFT_${pubField}-${tstart}-${Tsft}.sft"
    if ! test -f $MSFTsft; then
        echo "ERROR: could not find file '$MSFTsft'"
        exit 1
    fi

done
