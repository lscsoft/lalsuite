##---------- names of codes and input/output files
psd_code="lalapps_ComputePSD"
mfd_code="lalapps_Makefakedata_v4"

retstatus=0

# ---------- fixed parameters of our test SFTs
sqrtSh=3.25e-22
h0=5e-23
cosi=0
startTime=828002611
fmin=50
linefreq="50.05"
TSFT=1800

# ---------- fixed parameters for PSD runs
blocksRngMed=101

# ---------- define comparison functions
awk_absdiff='{printf "%.6f", sqrt(($1-$2)*($1-$2)) }'
awk_reldev='{printf "%.2e", sqrt(($1-$2)*($1-$2))/(0.5*($1+$2)) }'
awk_isgtr='{if($1>$2) {print "1"}}'
awk_nonpos='{if(0.0>=$1) {print "1"}}'

## ----- test correct cropping of normalized SFT power on a single SFT
echo "----------------------------------------------------------------------"
echo "STEP 1: testing --outputNormSFT frequency range..."
echo "----------------------------------------------------------------------"
echo

outSFT="./testpsd_sft_H1"

## ----- create SFTs
cmdline="${mfd_code} --IFO=H1 --outSingleSFT=1 --outSFTbname=$outSFT --startTime=$startTime --duration=$TSFT --fmin=$fmin --Band=0.1 --noiseSqrtSh=$sqrtSh --lineFeature=1 --h0=$h0 --cosi=$cosi --Freq=$linefreq --randSeed=1"

echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

outPSD_full="./psd_fullsft.dat"
outPSD_band="./psd_freqband.dat"

## ----- run computePSD once for full SFT
cmdline="${psd_code} --inputData=$outSFT --outputPSD=$outPSD_full --blocksRngMed=$blocksRngMed --outputNormSFT=1"

echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$psd_code' ..."
    exit 1
fi

## ----- run computePSD for a second time, with restricted band
cmdline="${psd_code} --inputData=$outSFT --outputPSD=$outPSD_band --blocksRngMed=$blocksRngMed --outputNormSFT=1  --Freq=50.03 --FreqBand=0.04"

echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$psd_code' ..."
    exit 1
fi

# get the result file row at maximum, which should recover the line feature
topline_psd_full=$(sort -nr -k2,2 $outPSD_full | head -1)
toppsd_full=$(echo $topline_psd_full | awk '{print $2}')
toppsdfreq_full=$(echo $topline_psd_full | awk '{print $1}')
topline_power_full=$(sort -nr -k3,3 $outPSD_full | head -1)
toppower_full=$(echo $topline_power_full | awk '{print $3}')
toppowerfreq_full=$(echo $topline_power_full | awk '{print $1}')

topline_psd_band=$(sort -nr -k2,2 $outPSD_band | head -1)
toppsd_band=$(echo $topline_psd_band | awk '{print $2}')
toppsdfreq_band=$(echo $topline_psd_band | awk '{print $1}')
topline_power_band=$(sort -nr -k3,3 $outPSD_band | head -1)
toppower_band=$(echo $topline_power_band | awk '{print $3}')
toppowerfreq_band=$(echo $topline_power_band | awk '{print $1}')

echo "Loudest bins:"
echo "==>  full SFT: PSD=$toppsd_full at $toppsdfreq_full Hz, normPower=$toppower_full at $toppowerfreq_full Hz"
echo "==>  freqband: PSD=$toppsd_band at $toppsdfreq_band Hz, normPower=$toppower_band at $toppowerfreq_band Hz"

tolerance_freq=0.0
tolerance_rel=0.001
freqdiff_full=$(echo $toppowerfreq_full $linefreq | awk "$awk_absdiff")
freqdiff_band=$(echo $toppowerfreq_band $linefreq | awk "$awk_absdiff")
freqdiff_runs=$(echo $toppowerfreq_full $toppowerfreq_band | awk "$awk_absdiff")
toppsd_reldev=$(echo $toppsd_full $toppsd_band | awk "$awk_reldev")
toppower_reldev=$(echo $toppower_full $toppower_band | awk "$awk_reldev")
echo "Injected line recovered in full run with offset $freqdiff_full Hz and in band run with $freqdiff_band Hz - difference between runs: $freqdiff_runs Hz."
echo "Relative difference in highest PSD value is $toppsd_reldev and in highest normPower value $toppower_reldev."

fail1=$(echo $freqdiff_full $tolerance_freq | awk "$awk_isgtr")
fail2=$(echo $freqdiff_band $tolerance_freq | awk "$awk_isgtr")
fail3=$(echo $freqdiff_runs $tolerance_freq | awk "$awk_isgtr")
fail4=$(echo $toppsd_reldev $tolerance_rel  | awk "$awk_isgtr")
fail5=$(echo $toppower_reldev $tolerance_rel  | awk "$awk_isgtr")
if [ "$fail1" -o "$fail2" -o "$fail3" -o "$fail4" -o "$fail5" ]; then
    echo " ==> FAILED"
    retstatus=1
else
    echo " ==> OK"
    echo
    echo "========== OK. Line feature successfully recovered. =========="
    echo
fi


echo "----------------------------------------------------------------------"
echo "STEP 2: testing math ops over SFTs and IFOs..."
echo "----------------------------------------------------------------------"
echo

tolerance_rel=0.000001

## ----- create SFTs
outSFTbname="./testpsd_multisfts"
IFOs=(H1 L1)
for X in 0 1; do
    outSFT="${outSFTbname}_${IFOs[$X]}"
    cmdline="${mfd_code} --IFO=${IFOs[$X]} --outSingleSFT=1 --outSFTbname=$outSFT --startTime=$startTime --duration=$(echo "2*$TSFT" | bc) --fmin=$fmin --Band=0.1 --noiseSqrtSh=$sqrtSh --lineFeature=1 --h0=$h0 --cosi=$cosi --Freq=$linefreq --randSeed=$X"
    echo $cmdline;
    if ! eval $cmdline; then
        echo "Error.. something failed when running '$mfd_code' ..."
        exit 1
    fi
done
echo

# function to de-spaghettify the loop over multiple mathops, IFOs and timestamps below
get_psd () { # expected argument order: mthop IFO startTime endTime extrArgs
    psdfile="psd_$2_$3-$4.dat"
    # here we only use the --IFO constraint for single IFOs;
    # currently the code would silently interpret H1L1 as H1 only!
    if [ ${#2} == 2 ]; then
        IFObit="--IFO=$2"
    else
        IFObit=""
    fi
    cmdline="${psd_code} ${IFObit} --inputData=${outSFTbname}* --outputPSD=$psdfile --blocksRngMed=$blocksRngMed --outputNormSFT=1 --startTime=$3 --endTime=$4 --PSDmthopSFTs=$1 --PSDmthopIFOs=$1 $5"
    echo $cmdline;
    if ! eval $cmdline; then
        echo "Error.. something failed when running '$psd_code' ..."
        exit 1
    fi
    firstline_psd=$(grep -v '^%' $psdfile | head -1)
    freq=$(echo $firstline_psd | awk '{print $1}')
    psd=$(echo $firstline_psd | awk '{print $2}')
    power=$(echo $firstline_psd | awk '{print $3}')
    firstline_psd=($freq $psd $power)
}

## ----- test the various mthops - for simplicity, use the same over SFTs and IFOs
# use the same ordering as the MathOpType enum:
mthops_awk=('{printf "%.6e", $1+$2}' # ARITHMETIC_SUM
            '{printf "%.6e", 0.5*($1+$2)}' # ARITHMETIC_MEAN
            '{printf "%.6e", 0.5*($1+$2)}' # ARITHMETIC_MEDIAN - note for 2 elements this falls back to mean!
            '{printf "%.6e", 1./(1./$1+1./$2)}' # HARMONIC_SUM
            '{printf "%.6e", 2./(1./$1+1./$2)}' # HARMONIC_MEAN
            '{printf "%.6e", 1./sqrt(1./($1*$1)+1./($2*$2))}' # POWERMINUS2_SUM
            '{printf "%.6e", 1./sqrt(0.5*(1./($1*$1)+1./($2*$2)))}' # POWERMINUS2_MEAN
            '{printf "%.6e", ($1<$2) ? $1 : $2}' # MINIMUM
            '{printf "%.6e", ($1>$2) ? $1 : $2}' # MAXIMUM
           )
for mthop in 0 1 2 3 4 5 6 7 8; do
    startTime2=$(echo "$startTime + $TSFT" | bc)
    endTime1=$(echo "$startTime + 1" | bc)
    endTime2=$(echo "$startTime2 + 1" | bc)
    echo "Running ComputePSD with PSDmthopSFTs=PSDmthopIFOs=$mthop..."
    for IFO in "H1" "L1" "H1L1"; do
        for minT in $startTime $startTime2; do
            for maxT in $endTime1 $endTime2; do
                if (( $maxT > $minT )); then
                    get_psd $mthop $IFO $minT $maxT
                    # storing the results in dynamically named variables
                    # so they can be accessed again outside the loop
                    psdvar="psd_${IFO}_${minT}_${maxT}"
                    powvar="power_${IFO}_${minT}_${maxT}"
                    declare $psdvar=${firstline_psd[1]}
                    declare $powvar=${firstline_psd[2]}
                fi
            done
        done
    done

    echo "--------- Compare results (PSDmthopSFTs=PSDmthopIFOs=$mthop) ----------------------------------------------------------------------------------------"
    echo "                	[T0,T1]		[T1,T2]		[T0,T2]		awk		(reldev) [Tolerance=$tolerance_rel]"
    failed=0
    # first check the over-SFTs mathop for each IFO (and multi-IFO)
    for IFO in "H1" "L1" "H1L1"; do
        psd1=$(eval "echo \$"$(echo "psd_${IFO}_${startTime}_${endTime1}")"")
        psd2=$(eval "echo \$"$(echo "psd_${IFO}_${startTime2}_${endTime2}")"")
        psd12=$(eval "echo \$"$(echo "psd_${IFO}_${startTime}_${endTime2}")"")
        psdawk=$(echo $psd1 $psd2 | awk "${mthops_awk[$mthop]}")
        psd_reldev=$(echo $psd12 $psdawk | awk "$awk_reldev")
        nonpos_psd1=$(echo $psd1 | awk "$awk_nonpos")
        nonpos_psd2=$(echo $psd2 | awk "$awk_nonpos")
        nonpos_psd12=$(echo $psd12 | awk "$awk_nonpos")
        failed_tolerance=$(echo $psd_reldev $tolerance_rel | awk "$awk_isgtr")
        if [ "$nonpos_psd1" -o "$nonpos_psd2" -o "$nonpos_psd12" -o "$failed_tolerance" ]; then
            failbit="FAILED!"
            failed=1
        else
            failbit="OK!"
        fi
        echo    "==>  $IFO:         	$psd1 	$psd2   	$psd12   	$psdawk	($psd_reldev) $failbit"
    done

    # then check the over-IFOs mathop (only for the full SFT set)
    psdH1=$(eval "echo \$"$(echo "psd_H1_${startTime}_${endTime2}")"")
    psdL1=$(eval "echo \$"$(echo "psd_L1_${startTime}_${endTime2}")"")
    psdH1L1=$(eval "echo \$"$(echo "psd_H1L1_${startTime}_${endTime2}")"")
    psdawk=$(echo $psdH1 $psdL1 | awk "${mthops_awk[$mthop]}")
    psd_reldev=$(echo $psdH1L1 $psdawk | awk "$awk_reldev")
    nonpos_psdH1L1=$(echo $psdH1L1 | awk "$awk_nonpos")
    failed_tolerance=$(echo $psd_reldev $tolerance_rel | awk "$awk_isgtr")
    if [ "$nonpos_psdH1L1" -o "$failed_tolerance" ]; then
        failbit="FAILED!"
        failed=1
    else
        failbit="OK!"
    fi
    echo    "==>  H1L1(awk):         	   	     	        $psdH1L1   	$psdawk	($psd_reldev) $failbit"

    if (( $failed > 0 )); then
        echo " ==> something FAILED"
        retstatus=1
    else
        echo "========== OK. All results consistent. =========="
        echo
    fi
done # loop over mthops


echo "----------------------------------------------------------------------"
echo "STEP 3: testing special case for combined all-IFO-all-SFTs harmmean/powerminus2mean..."
echo "----------------------------------------------------------------------"
echo


## ----- test the special normalizeByTotalNumSFTs mode
mthops_loop=('harmsum' # together with --normalizeByTotalNumSFTs should emulate harmmean
             'powerminus2sum' # together with --normalizeByTotalNumSFTs should emulate powerminus2mean
            )
mthops_awk=('{printf "%.6e", 4./(1./$1+1./$2+1./$3+1./$4)}' # harmonic mean over 4 inputs
            '{printf "%.6e", sqrt(4.)/sqrt((1./($1*$1)+1./($2*$2)+1./($3*$3)+1./($4*$4)))}' # power-minus-2 mean over 4 inputs
           )
for mthopidx in 0 1; do
    startTime2=$(echo "$startTime + $TSFT" | bc)
    endTime1=$(echo "$startTime + 1" | bc)
    endTime2=$(echo "$startTime2 + 1" | bc)
    echo "Running ComputePSD with PSDmthopSFTs=PSDmthopIFOs=${mthops_loop[$mthopidx]} and --normalizeByTotalNumSFTs..."
    for IFO in "H1" "L1" "H1L1"; do
        for minT in $startTime $startTime2; do
            for maxT in $endTime1 $endTime2; do
                if (( $maxT > $minT )); then
                    get_psd ${mthops_loop[$mthopidx]} $IFO $minT $maxT "--normalizeByTotalNumSFTs"
                    # storing the results in dynamically named variables
                    # so they can be accessed again outside the loop
                    psdvar="psd_${IFO}_${minT}_${maxT}"
                    powvar="power_${IFO}_${minT}_${maxT}"
                    declare $psdvar=${firstline_psd[1]}
                    declare $powvar=${firstline_psd[2]}
                fi
            done
        done
    done

    echo "--------- Compare results (PSDmthopSFTs=PSDmthopIFOs=${mthops_loop[$mthopidx]} --normalizeByTotalNumSFTs) ----------------------------------------------------------------------------------------"
    echo "                	[T0,T1]		[T1,T2]		[T0,T2]		awk		(reldev) [Tolerance=$tolerance_rel]"
    for IFO in "H1" "L1" "H1L1"; do
        psd1=$(eval "echo \$"$(echo "psd_${IFO}_${startTime}_${endTime1}")"")
        psd2=$(eval "echo \$"$(echo "psd_${IFO}_${startTime2}_${endTime2}")"")
        echo    "==>  $IFO:         	$psd1 	$psd2"
    done
    failed=0
    # check the combination over the full all-IFOs-all-SFTs set
    psdH11=$(eval "echo \$"$(echo "psd_H1_${startTime}_${endTime1}")"")
    psdH12=$(eval "echo \$"$(echo "psd_H1_${startTime2}_${endTime2}")"")
    psdL11=$(eval "echo \$"$(echo "psd_L1_${startTime}_${endTime1}")"")
    psdL12=$(eval "echo \$"$(echo "psd_L1_${startTime2}_${endTime2}")"")
    psdH1L1=$(eval "echo \$"$(echo "psd_H1L1_${startTime}_${endTime2}")"")
    psdawk=$(echo $psdH11 $psdH12 $psdL11 $psdL12 | awk "${mthops_awk[$mthopidx]}")
    psd_reldev=$(echo $psdH1L1 $psdawk | awk "$awk_reldev")
    nonpos_psdH1L1=$(echo $psdH1L1 | awk "$awk_nonpos")
    failed_tolerance=$(echo $psd_reldev $tolerance_rel | awk "$awk_isgtr")
    if [ "$nonpos_psdH1L1" -o "$failed_tolerance" ]; then
        failbit="FAILED!"
        failed=1
    else
        failbit="OK!"
    fi
    echo    "==>  H1L1(awk):         	   	     	        $psdH1L1   	$psdawk	($psd_reldev) $failbit"

    if (( $failed > 0 )); then
        echo " ==> something FAILED"
        retstatus=1
    else
        echo "========== OK. All results consistent. =========="
        echo
    fi
done # loop over mthops

exit $retstatus
