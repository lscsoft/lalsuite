if test "${LALFRAME_ENABLED}" = false; then
    echo "Skipping test: requires LALFrame"
    exit 77
fi

## common variables
tstart=1257741529
Tsft=1800
tend=$(echo "${tstart} + ${Tsft}" | bc)
Band=256

## run MFDv5 to create fake frame and equivalent SFT
mkdir -p MFDv5/
cmdline="lalpulsar_Makefakedata_v5 --IFOs H1 --sqrtSX 1e-24 --Tsft ${Tsft} --fmin 0 --Band ${Band} --injectionSources '{Alpha=0.1; Delta=0.4; Freq=50; f1dot=1e-10; h0=1e-24; cosi=0.7; refTime=${tstart}}' --outSingleSFT false --outFrameDir MFDv5/ --outSFTdir MFDv5/ --startTime ${tstart} --duration ${Tsft}"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
MFDv5gwf="./MFDv5/H-H1_mfdv5-${tstart}-${Tsft}.gwf"
if ! test -f $MFDv5gwf; then
    echo "ERROR: could not find file '$MFDv5gwf'"
    exit 1
fi
MFDv5sft="./MFDv5/H-1_H1_${Tsft}SFT_mfdv5-${tstart}-${Tsft}.sft"
if ! test -f $MFDv5sft; then
    echo "ERROR: could not find file '$MFDv5sft'"
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
mkdir -p MSFTs/
pubObsRun=4
pubObsKind="DEV"
pubRevision=1
for MCFch in AdcREAL4 AdcREAL8 ProcREAL4 ProcREAL8 SimREAL4 SimREAL8; do
    echo "=== frame channel: ${MCFch} ==="

    ## make a cache file for the MultiChFrame frames
    MCFframecache="./MCFframecache"
    echo "H H1_MultiChFrame ${tstart} 1800 file://localhost$PWD/MultiChFrame/H-H1_MultiChFrame-${tstart}-${Tsft}.gwf" > $MCFframecache

    ## run MakeSFTs to create SFTs from the fake frames
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

# check that SFTs made by MFDv5 vs from Adc channels (which are not scaled by makeMultiChFrame) are the same
tol=1e-10
for MCFch in AdcREAL4 AdcREAL8; do
    echo "=== frame channel: ${MCFch} ==="
    pubField="O${pubObsRun}${pubObsKind}+R${pubRevision}+C${MCFch}+WRECT"
    MSFTsft="./MSFTs/H-1_H1_${Tsft}SFT_${pubField}-${tstart}-${Tsft}.sft"
    cmdline="lalpulsar_compareSFTs -V -e $tol -1 $MFDv5sft -2 $MSFTsft"
    echo "Comparing SFTs produced by MFDv5 and MakeSFTs from ${MCFch} channel, allowed tolerance=$tol:"
    if ! eval "$cmdline"; then
        echo "ERROR: something failed when running '$cmdline'"
        exit 1
    fi
done

# check that SFTs made from REAL4 vs REAL8 channels (which are the same data) are also the same
for MCFchtype in Adc Proc Sim; do
    echo "=== frame channel type: ${MCFchtype} ==="
    pubFieldREAL4="O${pubObsRun}${pubObsKind}+R${pubRevision}+C${MCFchtype}REAL4+WRECT"
    MSFTsftREAL4="./MSFTs/H-1_H1_${Tsft}SFT_${pubFieldREAL4}-${tstart}-${Tsft}.sft"
    pubFieldREAL8="O${pubObsRun}${pubObsKind}+R${pubRevision}+C${MCFchtype}REAL8+WRECT"
    MSFTsftREAL8="./MSFTs/H-1_H1_${Tsft}SFT_${pubFieldREAL8}-${tstart}-${Tsft}.sft"
    cmdline="lalpulsar_compareSFTs -V -e $tol -1 $MSFTsftREAL4 -2 $MSFTsftREAL8"
    echo "Comparing SFTs produced by MakeSFTs from REAL4 vs REAL8 channels, allowed tolerance=$tol:"
    if ! eval "$cmdline"; then
        echo "ERROR: something failed when running '$cmdline'"
        exit 1
    fi
done

## run MakeSFTs to create multiple SFTs from different channels in the same frame
pubObsRun=4
pubObsKind="AUX"
pubRevision=1
cmdline="lalpulsar_MakeSFTs --frame-cache $MCFframecache --channel-name H1:AdcREAL4,H1:AdcREAL8,H1:ProcREAL4,H1:ProcREAL8,H1:SimREAL4,H1:SimREAL8 --sft-duration ${Tsft} --high-pass-freq 0 --start-freq 0 --band ${Band} --comment-field 'Test comment' --observing-run ${pubObsRun} --observing-kind ${pubObsKind} --observing-revision ${pubRevision} --sft-write-path MSFTs/,MSFTs/,MSFTs/,MSFTs/,MSFTs/,MSFTs/ --gps-start-time ${tstart} --gps-end-time ${tend} --window-type rectangular"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi

for MCFch in AdcREAL4 AdcREAL8 ProcREAL4 ProcREAL8 SimREAL4 SimREAL8; do

    echo "=== frame channel: ${MCFch} ==="
    pubFieldMulti="O${pubObsRun}${pubObsKind}+R${pubRevision}+C${MCFch}+WRECT"
    MSFTsftMulti="./MSFTs/H-1_H1_${Tsft}SFT_${pubFieldMulti}-${tstart}-${Tsft}.sft"
    if ! test -f $MSFTsftMulti; then
	echo "ERROR: could not find file '$MSFTsftMulti'"
	exit 1
    fi

    pubField="O${pubObsRun}DEV+R${pubRevision}+C${MCFch}+WRECT"
    MSFTsft="./MSFTs/H-1_H1_${Tsft}SFT_${pubField}-${tstart}-${Tsft}.sft"
    cmdline="lalpulsar_compareSFTs -V -e $tol -1 $MSFTsft -2 $MSFTsftMulti"
    echo "Comparing SFTs produced by MakeSFTs from single vs multiple channels, allowed tolerance=$tol:"
    if ! eval "$cmdline"; then
        echo "ERROR: something failed when running '$cmdline'"
        exit 1
    fi

done
