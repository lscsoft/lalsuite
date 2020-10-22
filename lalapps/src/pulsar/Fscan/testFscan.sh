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

## run spec_avg
outdir=spec_avg
mkdir -p ${outdir}
fmax=`echo ${fmin} + ${Band} | bc`
cmdline="( cd ${outdir} && ../lalapps_spec_avg -p ../${MFDv5sft} -I H1 -s 0 -e 2000000000 -f ${fmin} -F ${fmax} -t ${Tsft} -r 0.1 )"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
outfilebase="${outdir}/spec_10.00_20.00_H1_0_2000000000"
for file in ${outfilebase} ${outfilebase}_date ${outfilebase}_timeaverage ${outfilebase}_timestamps; do
    if ! test -f $file; then
        echo "ERROR: could not find file '$file'"
        exit 1
    fi
done

## run spec_avg_long to create an average ASD
outdir=spec_avg_long
mkdir -p ${outdir}
fmax=`echo ${fmin} + ${Band} | bc`
cmdline="( cd ${outdir} && ../lalapps_spec_avg_long -p ../${MFDv5sft} -I H1 -s 0 -e 2000000000 -f ${fmin} -F ${fmax} -t ${Tsft} )"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
outfilebase="${outdir}/spec_10.00_20.00_H1_0_2000000000"
for file in ${outfilebase}.txt; do
    if ! test -f $file; then
        echo "ERROR: could not find file '$file'"
        exit 1
    fi
done

## do a simple test
echo -n "Comparing first line of spec_avg_long/spec_10.00_20.00_H1_0_2000000000.txt to reference ... "
firstline=`awk 'NR == 1 { print }' spec_avg_long/spec_10.00_20.00_H1_0_2000000000.txt`
firstline_ref='10.00000000 2.1753e-45 4.66401e-23 1.84129e-45 4.29102e-23'
for f in 1 2 3 4 5; do
    field=`echo "$firstline" | awk '{ print $'"$f"' }'`
    field_ref=`echo "$firstline_ref" | awk '{ print $'"$f"' }'`
    cmdline="echo $field $field_ref | awk '{ exit ( sqrt( (\$1 - \$2)^2 ) / \$2 < 1e-5 ? 0 : 1 ) }'"
    if ! eval "$cmdline"; then
        echo "ERROR: $field and $field_ref differ by more than 1e-5 relative tolerance"
        exit 1
    fi
done
echo "OK"
echo -n "Comparing last line of spec_avg_long/spec_10.00_20.00_H1_0_2000000000.txt to reference ... "
lastline=`awk 'NR == 18001 { print }' spec_avg_long/spec_10.00_20.00_H1_0_2000000000.txt`
lastline_ref='20.00000000 6.68982e-45 8.17913e-23 6.43827e-45 8.02388e-23'
for f in 1 2 3 4 5; do
    field=`echo "$lastline" | awk '{ print $'"$f"' }'`
    field_ref=`echo "$lastline_ref" | awk '{ print $'"$f"' }'`
    cmdline="echo $field $field_ref | awk '{ exit ( sqrt( (\$1 - \$2)^2 ) / \$2 < 1e-5 ? 0 : 1 ) }'"
    if ! eval "$cmdline"; then
        echo "ERROR: $field and $field_ref differ by more than 1e-5 relative tolerance"
        exit 1
    fi
done
echo "OK"
