## common variables
tstart=1257741529
Tsft=1800
Band=256

## run MFDv5 to create a fake frame and equivalent SFT
cmdline="lalapps_Makefakedata_v5 --IFOs H1 --sqrtSX 1e-24 --startTime ${tstart} --duration ${Tsft} --fmin 0 --Band ${Band} --injectionSources '{Alpha=0.1; Delta=0.4; Freq=50; f1dot=1e-10; h0=1e-24; cosi=0.7; refTime=${tstart}}' --outFrameDir . --outSFTdir ."
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
MFDv5gwf="./H-H1_mfdv5-${tstart}-${Tsft}.gwf"
MFDv5sft="./H-1_H1_${Tsft}SFT_mfdv5-${tstart}-${Tsft}.sft"
for file in $MFDv5gwf $MFDv5sft; do
    if ! test -f $file; then
        echo "ERROR: could not find file '$file'"
        exit 1
    fi
done

## make a cache file for the MFDv5 frame
framecache="./framecache"
echo "H H1_mfdv5 ${tstart} ${Tsft} file://localhost$PWD/$MFDv5gwf" > $framecache

## run MakeSFTs to create an SFT from the fake frame
tend=`echo "${tstart} + ${Tsft}" | bc`
Band2=`echo "${Band} + 0.0006" | bc`   ## add 0.0006 to bandwidth to get same number of SFT bins as MDFv5
cmdline="lalapps_MakeSFTs -f 0 -t ${Tsft} -p . -C $framecache -s ${tstart} -e ${tend} -N H1:mfdv5 -v 2 -i H1 -u PROC_REAL8 -w 0 -F 0 -B ${Band2} -D 4 -X MSFT"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
MSFTsft="./H-1_H1_${Tsft}SFT_MSFT-1257/H-1_H1_${Tsft}SFT_MSFT-${tstart}-${Tsft}.sft"
for file in $MSFTsft; do
    if ! test -f $file; then
        echo "ERROR: could not find file '$file'"
        exit 1
    fi
done

## compare SFTs produced by MFDv5 and MakeSFTs, should be very nearly identical
tol=1e-10
cmdline="lalapps_compareSFTs -V -e $tol -1 $MFDv5sft -2 $MSFTsft"
echo "Comparing SFTs produced by MFDv5 and MakeSFTs, allowed tolerance=$tol:"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
