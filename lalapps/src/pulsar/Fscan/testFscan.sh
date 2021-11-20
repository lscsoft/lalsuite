## common variables
tstart=1000000000
dur=3600
Tsft=1800
fmin=10
Band=10
fmax=`echo ${fmin} + ${Band} | bc`

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

## do a simple test
echo -n "Comparing first line of spec_avg/spec_10.00_20.00_H1_0_2000000000_date to reference ... "
firstline=`awk 'NR == 1 { print }' spec_avg/spec_10.00_20.00_H1_0_2000000000_date`
firstline_ref='0	 2011	 9	 14	 1	 46	 25'
for f in 1 2 3 4 5 6 7; do
    field=`echo "$firstline" | awk '{ print $'"$f"' }'`
    field_ref=`echo "$firstline_ref" | awk '{ print $'"$f"' }'`
    cmdline="echo $field $field_ref | awk '{ exit ( sqrt( (\$1 - \$2)^2 ) < 1e-5 ? 0 : 1 ) }'"
    if ! eval "$cmdline"; then
        echo "ERROR: $field and $field_ref differ by more than 1e-5 relative tolerance"
        printf '=== firstline ===\n%s\n=== firstline_ref ===\n%s\n---' "${firstline}" "${firstline_ref}"
        exit 1
    fi
done
echo "OK"
echo -n "Comparing last line of spec_avg/spec_10.00_20.00_H1_0_2000000000_date to reference ... "
lastline=`awk 'NR == 2 { print }' spec_avg/spec_10.00_20.00_H1_0_2000000000_date`
lastline_ref='1	 2011	 9	 14	 2	 16	 25'
for f in 1 2 3 4 5 6 7; do
    field=`echo "$lastline" | awk '{ print $'"$f"' }'`
    field_ref=`echo "$lastline_ref" | awk '{ print $'"$f"' }'`
    cmdline="echo $field $field_ref | awk '{ exit ( sqrt( (\$1 - \$2)^2 ) < 1e-5 ? 0 : 1 ) }'"
    if ! eval "$cmdline"; then
        echo "ERROR: $field and $field_ref differ by more than 1e-5 relative tolerance"
        printf '=== lastline ===\n%s\n=== lastline_ref ===\n%s\n---' "${lastline}" "${lastline_ref}"
        exit 1
    fi
done
echo "OK"
echo -n "Comparing first line of spec_avg/spec_10.00_20.00_H1_0_2000000000_timestamps to reference ... "
firstline=`awk 'NR == 1 { print }' spec_avg/spec_10.00_20.00_H1_0_2000000000_timestamps`
firstline_ref='0.	1000000000'
for f in 1 2; do
    field=`echo "$firstline" | awk '{ print $'"$f"' }'`
    field_ref=`echo "$firstline_ref" | awk '{ print $'"$f"' }'`
    cmdline="echo $field $field_ref | awk '{ exit ( sqrt( (\$1 - \$2)^2 ) < 1e-5 ? 0 : 1 ) }'"
    if ! eval "$cmdline"; then
        echo "ERROR: $field and $field_ref differ by more than 1e-5 relative tolerance"
        printf '=== firstline ===\n%s\n=== firstline_ref ===\n%s\n---' "${firstline}" "${firstline_ref}"
        exit 1
    fi
done
echo "OK"
echo -n "Comparing last line of spec_avg/spec_10.00_20.00_H1_0_2000000000_timestamps to reference ... "
lastline=`awk 'NR == 2 { print }' spec_avg/spec_10.00_20.00_H1_0_2000000000_timestamps`
lastline_ref='1.	1000001800'
for f in 1 2; do
    field=`echo "$lastline" | awk '{ print $'"$f"' }'`
    field_ref=`echo "$lastline_ref" | awk '{ print $'"$f"' }'`
    cmdline="echo $field $field_ref | awk '{ exit ( sqrt( (\$1 - \$2)^2 ) < 1e-5 ? 0 : 1 ) }'"
    if ! eval "$cmdline"; then
        echo "ERROR: $field and $field_ref differ by more than 1e-5 relative tolerance"
        printf '=== lastline ===\n%s\n=== lastline_ref ===\n%s\n---' "${lastline}" "${lastline_ref}"
        exit 1
    fi
done
echo "OK"
echo -n "Comparing first line of spec_avg/spec_10.00_20.00_H1_0_2000000000_timeaverage to reference ... "
firstline=`awk 'NR == 1 { print }' spec_avg/spec_10.00_20.00_H1_0_2000000000_timeaverage`
firstline_ref='       10.000000            0.239'
for f in 1 2; do
    field=`echo "$firstline" | awk '{ print $'"$f"' }'`
    field_ref=`echo "$firstline_ref" | awk '{ print $'"$f"' }'`
    cmdline="echo $field $field_ref | awk '{ exit ( sqrt( (\$1 - \$2)^2 ) / \$2 < 1e-5 ? 0 : 1 ) }'"
    if ! eval "$cmdline"; then
        echo "ERROR: $field and $field_ref differ by more than 1e-5 relative tolerance"
        printf '=== firstline ===\n%s\n=== firstline_ref ===\n%s\n---' "${firstline}" "${firstline_ref}"
        exit 1
    fi
done
echo "OK"
echo -n "Comparing last line of spec_avg/spec_10.00_20.00_H1_0_2000000000_timeaverage to reference ... "
lastline=`awk 'NR == 18000 { print }' spec_avg/spec_10.00_20.00_H1_0_2000000000_timeaverage`
lastline_ref='       19.999444            1.224'
for f in 1 2; do
    field=`echo "$lastline" | awk '{ print $'"$f"' }'`
    field_ref=`echo "$lastline_ref" | awk '{ print $'"$f"' }'`
    cmdline="echo $field $field_ref | awk '{ exit ( sqrt( (\$1 - \$2)^2 ) / \$2 < 1e-5 ? 0 : 1 ) }'"
    if ! eval "$cmdline"; then
        echo "ERROR: $field and $field_ref differ by more than 1e-5 relative tolerance"
        printf '=== lastline ===\n%s\n=== lastline_ref ===\n%s\n---' "${lastline}" "${lastline_ref}"
        exit 1
    fi
done
echo "OK"
echo -n "Comparing first line of spec_avg/spec_10.00_20.00_H1_0_2000000000 to reference ... "
firstline=`awk 'NR == 1 { print }' spec_avg/spec_10.00_20.00_H1_0_2000000000`
firstline_ref='1.015133e-22    1.035036e-22    9.887402e-23    1.008830e-22'
for f in 1 2 3 4; do
    field=`echo "$firstline" | awk '{ print $'"$f"' }'`
    field_ref=`echo "$firstline_ref" | awk '{ print $'"$f"' }'`
    cmdline="echo $field $field_ref | awk '{ exit ( sqrt( (\$1 - \$2)^2 ) / \$2 < 1e-5 ? 0 : 1 ) }'"
    if ! eval "$cmdline"; then
        echo "ERROR: $field and $field_ref differ by more than 1e-5 relative tolerance"
        printf '=== firstline ===\n%s\n=== firstline_ref ===\n%s\n---' "${firstline}" "${firstline_ref}"
        exit 1
    fi
done
echo "OK"
echo -n "Comparing last line of spec_avg/spec_10.00_20.00_H1_0_2000000000 to reference ... "
lastline=`awk 'NR == 2 { print }' spec_avg/spec_10.00_20.00_H1_0_2000000000`
lastline_ref='9.351592e-23    9.947149e-23    9.312546e-23    9.680045e-23'
for f in 1 2 3 4; do
    field=`echo "$lastline" | awk '{ print $'"$f"' }'`
    field_ref=`echo "$lastline_ref" | awk '{ print $'"$f"' }'`
    cmdline="echo $field $field_ref | awk '{ exit ( sqrt( (\$1 - \$2)^2 ) / \$2 < 1e-5 ? 0 : 1 ) }'"
    if ! eval "$cmdline"; then
        echo "ERROR: $field and $field_ref differ by more than 1e-5 relative tolerance"
        printf '=== lastline ===\n%s\n=== lastline_ref ===\n%s\n---' "${lastline}" "${lastline_ref}"
        exit 1
    fi
done
echo "OK"

## run spec_avg_long to create an average ASD
outdir=spec_avg_long
mkdir -p ${outdir}
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
firstline_ref='10.00000000 2.1753e-45 4.66401e-23 1.84129e-45 4.29102e-23 0'
for f in 1 2 3 4 5 6; do
    field=`echo "$firstline" | awk '{ print $'"$f"' }'`
    field_ref=`echo "$firstline_ref" | awk '{ print $'"$f"' }'`
    if f < 6; then
       cmdline="echo $field $field_ref | awk '{ exit ( sqrt( (\$1 - \$2)^2 ) / \$2 < 1e-5 ? 0 : 1 ) }'"
    else
	cmdline="echo $field $field_ref | awk '{ exit ( sqrt( (\$1 - \$2)^2 ) < 1e-5 ? 0 : 1 ) }'"
    fi
    if ! eval "$cmdline"; then
        echo "ERROR: $field and $field_ref differ by more than 1e-5 relative tolerance"
        printf '=== firstline ===\n%s\n=== firstline_ref ===\n%s\n---' "${firstline}" "${firstline_ref}"
        exit 1
    fi
done
echo "OK"
echo -n "Comparing last line of spec_avg_long/spec_10.00_20.00_H1_0_2000000000.txt to reference ... "
lastline=`awk 'NR == 18000 { print }' spec_avg_long/spec_10.00_20.00_H1_0_2000000000.txt`
lastline_ref='19.99944444 1.26395e-44 1.12425e-22 1.29723e-44 1.13896e-22 0'
for f in 1 2 3 4 5 6; do
    field=`echo "$lastline" | awk '{ print $'"$f"' }'`
    field_ref=`echo "$lastline_ref" | awk '{ print $'"$f"' }'`
    if f < 6; then
       cmdline="echo $field $field_ref | awk '{ exit ( sqrt( (\$1 - \$2)^2 ) / \$2 < 1e-5 ? 0 : 1 ) }'"
    else
	cmdline="echo $field $field_ref | awk '{ exit ( sqrt( (\$1 - \$2)^2 ) < 1e-5 ? 0 : 1 ) }'"
    fi
    if ! eval "$cmdline"; then
        echo "ERROR: $field and $field_ref differ by more than 1e-5 relative tolerance"
        printf '=== lastline ===\n%s\n=== lastline_ref ===\n%s\n---' "${lastline}" "${lastline_ref}"
        exit 1
    fi
done
echo "OK"

## run spec_coherence to compute coherence values
outdir=spec_coherence
mkdir -p ${outdir}
cmdline="( cd ${outdir} && ../lalapps_spec_coherence -p ../${MFDv5sft} -q ../${MFDv5sft} -s 0 -e 2000000000 -f ${fmin} -F ${fmax} -t ${Tsft} )"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
outfilebase="${outdir}/spec_10.00_20.00_0_2000000000_coh"
for file in ${outfilebase}.txt; do
    if ! test -f $file; then
        echo "ERROR: could not find file '$file'"
        exit 1
    fi
done

## do a simple test
echo -n "Comparing first line of spec_coherence/spec_10.00_20.00_0_2000000000_coh.txt to reference ... "
firstline=`awk 'NR == 1 { print }' spec_coherence/spec_10.00_20.00_0_2000000000_coh.txt`
firstline_ref='10.00000000 1'
for f in 1 2; do
    field=`echo "$firstline" | awk '{ print $'"$f"' }'`
    field_ref=`echo "$firstline_ref" | awk '{ print $'"$f"' }'`
    cmdline="echo $field $field_ref | awk '{ exit ( sqrt( (\$1 - \$2)^2 ) / \$2 < 1e-5 ? 0 : 1 ) }'"
    if ! eval "$cmdline"; then
        echo "ERROR: $field and $field_ref differ by more than 1e-5 relative tolerance"
        printf '=== firstline ===\n%s\n=== firstline_ref ===\n%s\n---' "${firstline}" "${firstline_ref}"
        exit 1
    fi
done
echo "OK"
echo -n "Comparing last line of spec_coherence/spec_10.00_20.00_0_2000000000_coh.txt to reference ... "
lastline=`awk 'NR == 18000 { print }' spec_coherence/spec_10.00_20.00_0_2000000000_coh.txt`
lastline_ref='19.99944444 1'
for f in 1 2; do
    field=`echo "$lastline" | awk '{ print $'"$f"' }'`
    field_ref=`echo "$lastline_ref" | awk '{ print $'"$f"' }'`
    cmdline="echo $field $field_ref | awk '{ exit ( sqrt( (\$1 - \$2)^2 ) / \$2 < 1e-5 ? 0 : 1 ) }'"
    if ! eval "$cmdline"; then
        echo "ERROR: $field and $field_ref differ by more than 1e-5 relative tolerance"
        printf '=== lastline ===\n%s\n=== lastline_ref ===\n%s\n---' "${lastline}" "${lastline_ref}"
        exit 1
    fi
done
echo "OK"

## run MFDv5 to create some fake data
cmdline="lalapps_Makefakedata_v5 --outSingleSFT=TRUE --outSFTdir=. --IFOs=L1 --sqrtSX=1e-22 --startTime=${tstart} --duration=${dur} --fmin=${fmin} --Band=${Band} --Tsft=${Tsft} --SFToverlap=0 --randSeed=43"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
MFDv5sft2="./L-2_L1_${Tsft}SFT_mfdv5-${tstart}-${dur}.sft"
for file in $MFDv5sft2; do
    if ! test -f $file; then
        echo "ERROR: could not find file '$file'"
        exit 1
    fi
done

## run spec_coherence to compute coherence values
outdir=spec_coherence
mkdir -p ${outdir}
cmdline="( cd ${outdir} && ../lalapps_spec_coherence -p ../${MFDv5sft} -q ../${MFDv5sft2} -s 0 -e 2000000000 -f ${fmin} -F ${fmax} -t ${Tsft} )"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
outfilebase="${outdir}/spec_10.00_20.00_0_2000000000_coh"
for file in ${outfilebase}.txt; do
    if ! test -f $file; then
        echo "ERROR: could not find file '$file'"
        exit 1
    fi
done

## do a simple test
echo -n "Comparing first line of spec_coherence/spec_10.00_20.00_0_2000000000_coh.txt to reference ... "
firstline=`awk 'NR == 1 { print }' spec_coherence/spec_10.00_20.00_0_2000000000_coh.txt`
firstline_ref='10.00000000 0.0284384'
for f in 1 2; do
    field=`echo "$firstline" | awk '{ print $'"$f"' }'`
    field_ref=`echo "$firstline_ref" | awk '{ print $'"$f"' }'`
    cmdline="echo $field $field_ref | awk '{ exit ( sqrt( (\$1 - \$2)^2 ) / \$2 < 1e-5 ? 0 : 1 ) }'"
    if ! eval "$cmdline"; then
        echo "ERROR: $field and $field_ref differ by more than 1e-5 relative tolerance"
        printf '=== firstline ===\n%s\n=== firstline_ref ===\n%s\n---' "${firstline}" "${firstline_ref}"
        exit 1
    fi
done
echo "OK"
echo -n "Comparing last line of spec_coherence/spec_10.00_20.00_0_2000000000_coh.txt to reference ... "
lastline=`awk 'NR == 18000 { print }' spec_coherence/spec_10.00_20.00_0_2000000000_coh.txt`
lastline_ref='19.99944444 0.0888071'
for f in 1 2; do
    field=`echo "$lastline" | awk '{ print $'"$f"' }'`
    field_ref=`echo "$lastline_ref" | awk '{ print $'"$f"' }'`
    cmdline="echo $field $field_ref | awk '{ exit ( sqrt( (\$1 - \$2)^2 ) / \$2 < 1e-5 ? 0 : 1 ) }'"
    if ! eval "$cmdline"; then
        echo "ERROR: $field and $field_ref differ by more than 1e-5 relative tolerance"
        printf '=== lastline ===\n%s\n=== lastline_ref ===\n%s\n---' "${lastline}" "${lastline_ref}"
        exit 1
    fi
done
echo "OK"
