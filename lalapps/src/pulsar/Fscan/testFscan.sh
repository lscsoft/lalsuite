## common variables
tstart=1000000000
dur=3600
Tsft=1800
fmin=10
Band=10

## run MFDv5 to create some fake data
cmdline="lalapps_Makefakedata_v5 --outSingleSFT=TRUE --outSFTdir=. --IFOs=H1 --sqrtSX=1e-22 --startTime=${tstart} --duration=${dur} --fmin=${fmin} --Band=${Band} --Tsft=${Tsft} --SFToverlap=0 --randSeed=42"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
MFDv5sft="./H-2_H1_${Tsft}SFT_mfdv5-${tstart}-${dur}.sft"
for file in $MFDv5sft; do
    if ! test -f $file; then
        echo "ERROR: could not find file '$file'"
        exit 1
    fi
done

## run spec_avg_long to create an average ASD
fmax=`echo ${fmin} + ${Band} | bc`
cmdline="lalapps_spec_avg_long -p ${MFDv5sft} -I H1 -s 0 -e 2000000000 -f ${fmin} -F ${fmax} -t ${Tsft}"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
MSFTsft="./spec_10.00_20.00_H1_0_2000000000.txt"
for file in $MSFTsft; do
    if ! test -f $file; then
        echo "ERROR: could not find file '$file'"
        exit 1
    fi
done

## do a simple test
cmdline="head -n 1 spec_10.00_20.00_H1_0_2000000000.txt | awk '{if ($0=="     10.00000000 2.1753e-45 4.66401e-23 1.84129e-45 4.29102e-23") {print 1} else {print 0}}'"
echo "Comparing average spectra:"
if eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
