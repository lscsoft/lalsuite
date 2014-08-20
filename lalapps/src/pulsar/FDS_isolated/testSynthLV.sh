#!/bin/sh

## set LAL debug level
echo "Setting LAL_DEBUG_LEVEL=${LAL_DEBUG_LEVEL:-msglvl1,memdbg}"
export LAL_DEBUG_LEVEL

## allow 'make test' to work from builddir != srcdir
if [ -z "${srcdir}" ]; then
    srcdir=`dirname $0`
fi

## make sure we work in 'C' locale here to avoid awk sillyness
LC_ALL_old=$LC_ALL
export LC_ALL=C

builddir="./";
injectdir="../Injections/"
capdir="../Tools/"

##---------- names of codes
synth_code="${builddir}lalapps_synthesizeLVStats"
cap_code="${capdir}lalapps_ComputeAntennaPattern"
mfd_code="${injectdir}lalapps_Makefakedata_v4"
pfs_code="${builddir}lalapps_PredictFStat"

testDir="./testSynthLV_dir";
if [ -d "$testDir" ]; then
    rm -rf $testDir
fi
mkdir -p "$testDir"

## awk commands needed for testing
awk_absdiff='{printf "%.6f", sqrt(($1-$2)*($1-$2)) }'
awk_reldev='{printf "%.2e", sqrt(($1-$2)*($1-$2))/(0.5*sqrt(($1+$2)*($1+$2))) }'
awk_isgtr='{if($1>$2) {print "1"}}'
## awk line filtering: !/%%/ filters out header lines and /[0-9]/ filters out blank lines (require at least one number)
awk_print_wo_headers='!/%%/ && /[0-9]/ {print $col}'
awk_avg='!/%%/ && /[0-9]/ { total += $col; count++ } END { print total/count }'

# ---------- common test parameters
timestamp1=818845553
duration=90000
Tsft=1800
Nsteps=$(echo $duration $Tsft | awk '{print $1/$2}')
SNR=2
log10e=0.4342944819032518276511289189166051

## convergence is very slow, so optionally do many more draws and a bit tighter tolerance
if [ -n "$LOWTOL" ]; then
 numDraws=100000
 tolerance_F=0.025
 tolerance_p=0.025
 tolerance_LR=1e-5
else
 numDraws=1000
 tolerance_F=0.05
 tolerance_p=0.1
 tolerance_LR=1e-5
fi

outCAP1="${testDir}/antenna_pattern_test.dat"
outCAP2="${testDir}/antenna_pattern_test_n.dat"
stats_file_g="${testDir}/synthLV_test_stats_g.dat"
params_file_g="${testDir}/synthLV_test_params_g.dat"
stats_file_s="${testDir}/synthLV_test_stats_s.dat"
params_file_s="${testDir}/synthLV_test_params_s.dat"
stats_file_l="${testDir}/synthLV_test_stats_l.dat"
params_file_l="${testDir}/synthLV_test_params_l.dat"
stats_file_ns="${testDir}/synthLV_test_stats_ns.dat"
params_file_ns="${testDir}/synthLV_test_params_ns.dat"
stats_file_nlh="${testDir}/synthLV_test_stats_nlh.dat"
params_file_nlh="${testDir}/synthLV_test_params_nlh.dat"
stats_file_nll="${testDir}/synthLV_test_stats_nll.dat"
params_file_nll="${testDir}/synthLV_test_params_nll.dat"
timestampsfile="${testDir}/timestamps_test.dat"

# do all-sky average by default to get expected SNR scaling, use defaults for Fstar0 and oLGX for now
synth_cmdline_common="${synth_code} --IFOs=H1,L1 --outputMmunuX --dataStartGPS=$timestamp1 --dataDuration=$duration --numDraws=$numDraws --randSeed=1 --computeLR"

echo "----------------------------------------------------------------------------------------------------"
echo "synthesizeLVstats Test1a: average F-stats for pure Gaussian noise, signals, lines";
echo "----------------------------------------------------------------------------------------------------"

synth_cmdline="${synth_cmdline_common} --fixedSNR=0 --outputStats=$stats_file_g --outputInjParams=$params_file_g"
echo $synth_cmdline;
if ! eval $synth_cmdline; then
    echo "Error.. something failed when running '$synth_cmdline' ..."
    exit 1
fi

exp_2Fg=4.00

avg_2Fg=$(awk -v col=5 "$awk_avg" $stats_file_g)
reldev_avg_2Fg=$(echo $avg_2Fg $exp_2Fg | awk "$awk_reldev")
fail_avg_2Fg=$(echo $reldev_avg_2Fg $tolerance_F | awk "$awk_isgtr")

avg_2FH1g=$(awk -v col=6 "$awk_avg" $stats_file_g)
reldev_avg_2FH1g=$(echo $avg_2FH1g $exp_2Fg | awk "$awk_reldev")
fail_avg_2FH1g=$(echo $reldev_avg_2FH1g $tolerance_F | awk "$awk_isgtr")

avg_2FL1g=$(awk -v col=7 "$awk_avg" $stats_file_g)
reldev_avg_2FL1g=$(echo $avg_2FL1g $exp_2Fg | awk "$awk_reldev")
fail_avg_2FL1g=$(echo $reldev_avg_2FL1g $tolerance_F | awk "$awk_isgtr")

echo "==> average stats from $numDraws pure Gaussian noise draws:"
echo "    <2F>   =$avg_2Fg (expected: $exp_2Fg, rel. deviation: $reldev_avg_2Fg)"
echo "    <2F_H1>=$avg_2FH1g (expected: $exp_2Fg, rel. deviation: $reldev_avg_2FH1g)"
echo "    <2F_L1>=$avg_2FL1g (expected: $exp_2Fg, rel. deviation: $reldev_avg_2FL1g)"

if [ "$fail_avg_2Fg" -o "$fail_avg_2FH1g" -o "$fail_avg_2FL1g" ]; then
    echo "==> FAILED at tolerance=$tolerance_F"
    exit 1
else
    echo "==> OK at tolerance=$tolerance_F"
fi

synth_cmdline="${synth_cmdline_common} --fixedSNR=$SNR --outputStats=$stats_file_s --outputInjParams=$params_file_s"
echo $synth_cmdline;
if ! eval $synth_cmdline; then
    echo "Error.. something failed when running '$synth_cmdline' ..."
    exit 1
fi

exp_2Fs=$(echo $SNR | awk '{printf "%.2f",4+$1*$1}')
exp_2FH1s=$(echo $SNR | awk '{printf "%.2f",4+$1*$1/2}')
exp_2FL1s=$exp_2FH1s

avg_2Fs=$(awk -v col=5 "$awk_avg" $stats_file_s)
reldev_avg_2Fs=$(echo $avg_2Fs $exp_2Fs | awk "$awk_reldev")
fail_avg_2Fs=$(echo $reldev_avg_2Fs $tolerance_F | awk "$awk_isgtr")

avg_2FH1s=$(awk -v col=6 "$awk_avg" $stats_file_s)
reldev_avg_2FH1s=$(echo $avg_2FH1s $exp_2FH1s | awk "$awk_reldev")
fail_avg_2FH1s=$(echo $reldev_avg_2FH1s $tolerance_F | awk "$awk_isgtr")

avg_2FL1s=$(awk -v col=7 "$awk_avg" $stats_file_s)
reldev_avg_2FL1s=$(echo $avg_2FL1s $exp_2FL1s | awk "$awk_reldev")
fail_avg_2FL1s=$(echo $reldev_avg_2FL1s $tolerance_F | awk "$awk_isgtr")

echo "==> average stats from $numDraws draws with a SNR=2 signal:"
echo "    <2F>   =$avg_2Fs (expected: $exp_2Fs, rel. deviation: $reldev_avg_2Fs)"
echo "    <2F_H1>=$avg_2FH1s (expected: $exp_2FH1s, rel. deviation: $reldev_avg_2FH1s)"
echo "    <2F_L1>=$avg_2FL1s (expected: $exp_2FL1s, rel. deviation: $reldev_avg_2FL1s)"

if [ "$fail_avg_2Fs" -o "$fail_avg_2FH1s" -o "$fail_avg_2FL1s" ]; then
    echo "==> FAILED at tolerance=$tolerance_F"
    exit 1
else
    echo "==> OK at tolerance=$tolerance_F"
fi

synth_cmdline="${synth_cmdline_common} --fixedSNR=$SNR --lineIFO=L1 --outputStats=$stats_file_l --outputInjParams=$params_file_l"
echo $synth_cmdline;
if ! eval $synth_cmdline; then
    echo "Error.. something failed when running '$synth_cmdline' ..."
    exit 1
fi

exp_2Fl=$(echo $SNR | awk '{printf "%.2f",4+$1*$1/2}')
exp_2FH1l=4.00
exp_2FL1l=$(echo $SNR | awk '{printf "%.2f",4+$1*$1}')

avg_2Fl=$(awk -v col=5 "$awk_avg" $stats_file_l)
reldev_avg_2Fl=$(echo $avg_2Fl $exp_2Fl | awk "$awk_reldev")
fail_avg_2Fl=$(echo $reldev_avg_2Fl $tolerance_F | awk "$awk_isgtr")

avg_2FH1l=$(awk -v col=6 "$awk_avg" $stats_file_l)
reldev_avg_2FH1l=$(echo $avg_2FH1l $exp_2FH1l | awk "$awk_reldev")
fail_avg_2FH1l=$(echo $reldev_avg_2FH1l $tolerance_F | awk "$awk_isgtr")

avg_2FL1l=$(awk -v col=7 "$awk_avg" $stats_file_l)
reldev_avg_2FL1l=$(echo $avg_2FL1l $exp_2FL1l | awk "$awk_reldev")
fail_avg_2FL1l=$(echo $reldev_avg_2FL1l $tolerance_F | awk "$awk_isgtr")

echo "==> average stats from $numDraws draws with a SNR=2 line in L1:"
echo "    <2F>   =$avg_2Fl (expected: $exp_2Fl, rel. deviation: $reldev_avg_2Fl)"
echo "    <2F_H1>=$avg_2FH1l (expected: $exp_2FH1l, rel. deviation: $reldev_avg_2FH1l)"
echo "    <2F_L1>=$avg_2FL1l (expected: $exp_2FL1l, rel. deviation: $reldev_avg_2FL1l)"

if [ "$fail_avg_2Fl" -o "$fail_avg_2FH1l" -o "$fail_avg_2FL1l" ]; then
    echo "==> FAILED at tolerance=$tolerance_F"
    exit 1
else
    echo "==> OK at tolerance=$tolerance_F"
fi

echo "----------------------------------------------------------------------------------------------------"
echo "synthesizeLVstats Test1b: parameter consistency from previous draws";
echo "----------------------------------------------------------------------------------------------------"

## TODO: add a test here checking that params for g s l are all identical, which is assumed in the following

pi=3.1415926536
exp_alpha=$pi
exp_delta=0
exp_cosi=0
exp_psi=0
exp_phi0=$pi
avg_alpha=$(awk -v col=1 "$awk_avg" $params_file_g)
avg_delta=$(awk -v col=2 "$awk_avg" $params_file_g)
avg_cosi=$(awk -v col=5 "$awk_avg" $params_file_g)
avg_psi=$(awk -v col=6 "$awk_avg" $params_file_g)
avg_phi0=$(awk -v col=7 "$awk_avg" $params_file_g)
absdiff_alpha=$(echo $avg_alpha $exp_alpha | awk "$awk_absdiff")
absdiff_delta=$(echo $avg_delta $exp_delta | awk "$awk_absdiff")
absdiff_cosi=$(echo $avg_cosi $exp_cosi | awk "$awk_absdiff")
absdiff_psi=$(echo $avg_psi $exp_psi | awk "$awk_absdiff")
absdiff_phi0=$(echo $avg_phi0 $exp_phi0 | awk "$awk_absdiff")
fail_alpha=$(echo $absdiff_alpha $tolerance_p | awk "$awk_isgtr")
fail_delta=$(echo $absdiff_delta $tolerance_p | awk "$awk_isgtr")
fail_cosi=$(echo $absdiff_cosi $tolerance_p | awk "$awk_isgtr")
fail_psi=$(echo $absdiff_psi $tolerance_p | awk "$awk_isgtr")
fail_phi0=$(echo $absdiff_phi0 $tolerance_p | awk "$awk_isgtr")

echo "Parameter averages from previous draws: <alpha>=$avg_alpha, <delta>=$avg_delta, <cosi>=$avg_cosi, <psi>=$avg_psi, <phi0>=$avg_phi0"
echo "(expected:                                      $exp_alpha, $exp_delta, $exp_cosi, $exp_psi, $exp_phi0)"
echo "(abs.diff.:                                     $absdiff_alpha, $absdiff_delta, $absdiff_cosi, $absdiff_psi, $absdiff_phi0)"

if [ "$fail_alpha" -o "$fail_delta" -o "$fail_cosi" -o "$fail_psi" -o "$fail_phi0" ]; then
    echo "==> FAILED at tolerance=$tolerance_p"
    exit 1
else
    echo "==> OK at tolerance=$tolerance_p"
fi

echo "----------------------------------------------------------------------------------------------------"
echo "synthesizeLVstats Test1c: antenna patterns";
echo "----------------------------------------------------------------------------------------------------"

## ----- make timestampsfile
if [ -f $timestampsfile ]; then
    rm $timestampsfile
fi
iTS=1
while [ $iTS -le $Nsteps ]; do
    timestamp_i=$(echo $timestamp1 $iTS $Tsft| awk '{print $1 + ($2 - 1) * $3}')
    printf "%s 0\n" "$timestamp_i" >> $timestampsfile
    iTS=$(($iTS + 1))
done

params_draw1=$(awk '!/%%/ && /[0-9]/' $params_file_s | head -n 1)
alpha1=$(echo $params_draw1 | awk '{print $1}')
delta1=$(echo $params_draw1 | awk '{print $2}')
A_synth=$(echo $params_draw1 $Nsteps | awk '{print $12/(2*$28)}')
B_synth=$(echo $params_draw1 $Nsteps | awk '{print $13/(2*$28)}')
C_synth=$(echo $params_draw1 $Nsteps | awk '{print $14/(2*$28)}')
D_synth=$(echo $A_synth $B_synth $C_synth | awk '{print $1*$2-$3*$3}')
A_H1_synth=$(echo $params_draw1 $Nsteps | awk '{print $20/$28}')
B_H1_synth=$(echo $params_draw1 $Nsteps | awk '{print $21/$28}')
C_H1_synth=$(echo $params_draw1 $Nsteps | awk '{print $22/$28}')
D_H1_synth=$(echo $params_draw1 $Nsteps | awk '{print $23/($28*$28)}')
A_L1_synth=$(echo $params_draw1 $Nsteps | awk '{print $24/$28}')
B_L1_synth=$(echo $params_draw1 $Nsteps | awk '{print $25/$28}')
C_L1_synth=$(echo $params_draw1 $Nsteps | awk '{print $26/$28}')
D_L1_synth=$(echo $params_draw1 $Nsteps | awk '{print $27/($28*$28)}')

## ----- run ComputeAntennaPattern
cap_cmdline="${cap_code} --IFOs=H1,L1 --timeStampsFile=$timestampsfile --outABCD=$outCAP1 --Alpha=$alpha1 --Delta=$delta1"
echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi

A_cap=$(awk -v col=3 "$awk_print_wo_headers" $outCAP1)
B_cap=$(awk -v col=4 "$awk_print_wo_headers" $outCAP1)
C_cap=$(awk -v col=5 "$awk_print_wo_headers" $outCAP1)
D_cap=$(awk -v col=6 "$awk_print_wo_headers" $outCAP1)
A_H1_cap=$(awk -v col=7 "$awk_print_wo_headers" $outCAP1)
B_H1_cap=$(awk -v col=8 "$awk_print_wo_headers" $outCAP1)
C_H1_cap=$(awk -v col=9 "$awk_print_wo_headers" $outCAP1)
D_H1_cap=$(awk -v col=10 "$awk_print_wo_headers" $outCAP1)
A_L1_cap=$(awk -v col=11 "$awk_print_wo_headers" $outCAP1)
B_L1_cap=$(awk -v col=12 "$awk_print_wo_headers" $outCAP1)
C_L1_cap=$(awk -v col=13 "$awk_print_wo_headers" $outCAP1)
D_L1_cap=$(awk -v col=14 "$awk_print_wo_headers" $outCAP1)
reldev_A=$(echo $A_cap $A_synth | awk "$awk_reldev")
reldev_B=$(echo $B_cap $B_synth | awk "$awk_reldev")
reldev_C=$(echo $C_cap $C_synth | awk "$awk_reldev")
reldev_D=$(echo $D_cap $D_synth | awk "$awk_reldev")
reldev_A_H1=$(echo $A_H1_cap $A_H1_synth | awk "$awk_reldev")
reldev_B_H1=$(echo $B_H1_cap $B_H1_synth | awk "$awk_reldev")
reldev_C_H1=$(echo $C_H1_cap $C_H1_synth | awk "$awk_reldev")
reldev_D_H1=$(echo $D_H1_cap $D_H1_synth | awk "$awk_reldev")
reldev_A_L1=$(echo $A_L1_cap $A_L1_synth | awk "$awk_reldev")
reldev_B_L1=$(echo $B_L1_cap $B_L1_synth | awk "$awk_reldev")
reldev_C_L1=$(echo $C_L1_cap $C_L1_synth | awk "$awk_reldev")
reldev_D_L1=$(echo $D_L1_cap $D_L1_synth | awk "$awk_reldev")
fail_A=$(echo $reldev_A $tolerance_p | awk "$awk_isgtr")
fail_B=$(echo $reldev_B $tolerance_p | awk "$awk_isgtr")
fail_C=$(echo $reldev_C $tolerance_p | awk "$awk_isgtr")
fail_D=$(echo $reldev_D $tolerance_p | awk "$awk_isgtr")
fail_A_H1=$(echo $reldev_A_H1 $tolerance_p | awk "$awk_isgtr")
fail_B_H1=$(echo $reldev_B_H1 $tolerance_p | awk "$awk_isgtr")
fail_C_H1=$(echo $reldev_C_H1 $tolerance_p | awk "$awk_isgtr")
fail_D_H1=$(echo $reldev_D_H1 $tolerance_p | awk "$awk_isgtr")
fail_A_L1=$(echo $reldev_A_L1 $tolerance_p | awk "$awk_isgtr")
fail_B_L1=$(echo $reldev_B_L1 $tolerance_p | awk "$awk_isgtr")
fail_C_L1=$(echo $reldev_C_L1 $tolerance_p | awk "$awk_isgtr")
fail_D_L1=$(echo $reldev_D_L1 $tolerance_p | awk "$awk_isgtr")

echo "==> from lalapps_ComputeAntennaPattern: <A>=$A_cap, <B>=$B_cap, <C>=$C_cap, <D>=$D_cap, <A_H1>=$A_H1_cap, <B_H1>=$B_H1_cap, <C_H1>=$C_H1_cap, <D_H1>=$D_H1_cap, <A_L1>=$A_L1_cap, <B_L1>=$B_L1_cap, <C_L1>=$C_L1_cap, <D_L1>=$D_L1_cap"
echo "    from lalapps_synthesizeLVStats:     <A>=$A_synth, <B>=$B_synth, <C>=$C_synth, <D>=$D_synth, <A_H1>=$A_H1_synth, <B_H1>=$B_H1_synth, <C_H1>=$C_H1_synth, <D_H1>=$D_H1_synth, <A_L1>=$A_L1_synth, <B_L1>=$B_L1_synth, <C_L1>=$C_L1_synth, <D_L1>=$D_L1_synth"

if [ "$fail_A" -o "$fail_B" -o "$fail_C" -o "$fail_D" -o "$fail_A_H1" -o "$fail_B_H1" -o "$fail_C_H1" -o "$fail_D_H1" -o "$fail_A_L1" -o "$fail_B_L1" -o "$fail_C_L1" -o "$fail_D_L1" ]; then
    echo "==> FAILED at tolerance=$tolerance_p"
    exit 1
else
    echo "==> OK at tolerance=$tolerance_p"
fi

echo "----------------------------------------------------------------------------------------------------"
echo "synthesizeLVstats Test1d: LR-stats";
echo "----------------------------------------------------------------------------------------------------"

stats_draw1=$(awk '!/%%/ && /[0-9]/' $stats_file_s | head -n 1)
draw1_2F=$(echo $stats_draw1 | awk '{print $5}')
draw1_2FH1=$(echo $stats_draw1 | awk '{print $6}')
draw1_2FL1=$(echo $stats_draw1 | awk '{print $7}')
draw1_LR=$(echo $stats_draw1 | awk '{print $8}')

draw1_LR_recomp=$(echo $draw1_2F $draw1_2FH1 $draw1_2FL1 $log10e | awk '{print ( 0.5*$1 - log( 0.5 * ( exp(0.5*$2)*0.5 + exp(0.5*$3)*0.5 ) ) ) * $4}')
reldev_LR_draw1=$(echo $draw1_LR $draw1_LR_recomp | awk "$awk_reldev")
fail_LR_draw1=$(echo $reldev_LR_draw1 $tolerance_LR | awk "$awk_isgtr")

echo "for draw1, 2F=$draw1_2F, 2FH1=$draw1_2FH1, 2FL1=$draw1_2FL1:"
echo "           LR    =$draw1_LR  -- recomputed externally: $draw1_LR_recomp -- rel. deviation: $reldev_LR_draw1)"

if [ "$fail_LR_draw1" ]; then
    echo "==> FAILED at tolerance=$tolerance_LR"
    exit 1
else
    echo "==> OK at tolerance=$tolerance_LR"
fi


echo "----------------------------------------------------------------------------------------------------"
echo "synthesizeLVstats Test2a: F-stats for varying detector sensitivity";
echo "----------------------------------------------------------------------------------------------------"

Fstar0=10
oLGH1=0.2 # use weird values here to test internal conversions
oLGL1=3
oLG=$(echo $oLGH1 $oLGL1 | awk '{print $1+$2}')
pL=$(echo $oLG | awk '{print $1/(1+$1)}')
rH1=$(echo $oLGH1 $oLG | awk '{print $1/$2}') # these are not really the rX of Keitel et al 2014, but divided by Ndet
rL1=$(echo $oLGL1 $oLG | awk '{print $1/$2}')

synth_cmdline_common2="${synth_cmdline_common} --sqrtSX=1,2 --LR_Fstar0=$Fstar0  --LR_oLGX=$oLGH1,$oLGL1"

synth_cmdline="${synth_cmdline_common2} --fixedSNR=$SNR --outputStats=$stats_file_ns --outputInjParams=$params_file_ns"
echo $synth_cmdline;
if ! eval $synth_cmdline; then
    echo "Error.. something failed when running '$synth_cmdline' ..."
    exit 1
fi

exp_2Fs=$(echo $SNR | awk '{printf "%.2f",4+$1*$1}')
exp_2FH1s=$(echo $SNR | awk '{printf "%.2f",4+$1*$1*4/5}')
exp_2FL1s=$(echo $SNR | awk '{printf "%.2f",4+$1*$1/5}')

avg_2Fs=$(awk -v col=5 "$awk_avg" $stats_file_ns)
reldev_avg_2Fs=$(echo $avg_2Fs $exp_2Fs | awk "$awk_reldev")
fail_avg_2Fs=$(echo $reldev_avg_2Fs $tolerance_F | awk "$awk_isgtr")

avg_2FH1s=$(awk -v col=6 "$awk_avg" $stats_file_ns)
reldev_avg_2FH1s=$(echo $avg_2FH1s $exp_2FH1s | awk "$awk_reldev")
fail_avg_2FH1s=$(echo $reldev_avg_2FH1s $tolerance_F | awk "$awk_isgtr")

avg_2FL1s=$(awk -v col=7 "$awk_avg" $stats_file_ns)
reldev_avg_2FL1s=$(echo $avg_2FL1s $exp_2FL1s | awk "$awk_reldev")
fail_avg_2FL1s=$(echo $reldev_avg_2FL1s $tolerance_F | awk "$awk_isgtr")

echo "==> average stats from $numDraws draws with a SNR=2 signal (H1 2 times more sensitive than L1):"
echo "    <2F>   =$avg_2Fs (expected: $exp_2Fs, rel. deviation: $reldev_avg_2Fs)"
echo "    <2F_H1>=$avg_2FH1s (expected: $exp_2FH1s, rel. deviation: $reldev_avg_2FH1s)"
echo "    <2F_L1>=$avg_2FL1s (expected: $exp_2FL1s, rel. deviation: $reldev_avg_2FL1s)"

if [ "$fail_avg_2Fs" -o "$fail_avg_2FH1s" -o "$fail_avg_2FL1s" ]; then
    echo "==> FAILED at tolerance=$tolerance_F"
    exit 1
else
    echo "==> OK at tolerance=$tolerance_F"
fi

synth_cmdline="${synth_cmdline_common2} --fixedSNR=$SNR --outputStats=$stats_file_nlh --outputInjParams=$params_file_nlh --lineIFO=H1"
echo $synth_cmdline;
if ! eval $synth_cmdline; then
    echo "Error.. something failed when running '$synth_cmdline' ..."
    exit 1
fi

exp_2Fl=$(echo $SNR | awk '{printf "%.2f",4+$1*$1*4/5}')
exp_2FH1l=$(echo $SNR | awk '{printf "%.2f",4+$1*$1}')
exp_2FL1l=4.00

avg_2Fl=$(awk -v col=5 "$awk_avg" $stats_file_nlh)
reldev_avg_2Fl=$(echo $avg_2Fl $exp_2Fl | awk "$awk_reldev")
fail_avg_2Fl=$(echo $reldev_avg_2Fl $tolerance_F | awk "$awk_isgtr")

avg_2FH1l=$(awk -v col=6 "$awk_avg" $stats_file_nlh)
reldev_avg_2FH1l=$(echo $avg_2FH1l $exp_2FH1l | awk "$awk_reldev")
fail_avg_2FH1l=$(echo $reldev_avg_2FH1l $tolerance_F | awk "$awk_isgtr")

avg_2FL1l=$(awk -v col=7 "$awk_avg" $stats_file_nlh)
reldev_avg_2FL1l=$(echo $avg_2FL1l $exp_2FL1l | awk "$awk_reldev")
fail_avg_2FL1l=$(echo $reldev_avg_2FL1l $tolerance_F | awk "$awk_isgtr")

echo "==> average stats from $numDraws draws with a SNR=2 line in H1 (H1 2 times more sensitive than L1):"
echo "    <2F>   =$avg_2Fl (expected: $exp_2Fl, rel. deviation: $reldev_avg_2Fl)"
echo "    <2F_H1>=$avg_2FH1l (expected: $exp_2FH1l, rel. deviation: $reldev_avg_2FH1l)"
echo "    <2F_L1>=$avg_2FL1l (expected: $exp_2FL1l, rel. deviation: $reldev_avg_2FL1l)"

if [ "$fail_avg_2Fl" -o "$fail_avg_2FH1l" -o "$fail_avg_2FL1l" ]; then
    echo "==> FAILED at tolerance=$tolerance_F"
    exit 1
else
    echo "==> OK at tolerance=$tolerance_F"
fi

synth_cmdline="${synth_cmdline_common2} --fixedSNR=$SNR --outputStats=$stats_file_nll --outputInjParams=$params_file_nll --lineIFO=L1"
echo $synth_cmdline;
if ! eval $synth_cmdline; then
    echo "Error.. something failed when running '$synth_cmdline' ..."
    exit 1
fi

exp_2Fl=$(echo $SNR | awk '{printf "%.2f",4+$1*$1/5}')
exp_2FH1l=4.00
exp_2FL1l=$(echo $SNR | awk '{printf "%.2f",4+$1*$1}')

avg_2Fl=$(awk -v col=5 "$awk_avg" $stats_file_nll)
reldev_avg_2Fl=$(echo $avg_2Fl $exp_2Fl | awk "$awk_reldev")
fail_avg_2Fl=$(echo $reldev_avg_2Fl $tolerance_F | awk "$awk_isgtr")

avg_2FH1l=$(awk -v col=6 "$awk_avg" $stats_file_nll)
reldev_avg_2FH1l=$(echo $avg_2FH1l $exp_2FH1l | awk "$awk_reldev")
fail_avg_2FH1l=$(echo $reldev_avg_2FH1l $tolerance_F | awk "$awk_isgtr")

avg_2FL1l=$(awk -v col=7 "$awk_avg" $stats_file_nll)
reldev_avg_2FL1l=$(echo $avg_2FL1l $exp_2FL1l | awk "$awk_reldev")
fail_avg_2FL1l=$(echo $reldev_avg_2FL1l $tolerance_F | awk "$awk_isgtr")

echo "==> average stats from $numDraws draws with a SNR=2 line in L1 (H1 2 times more sensitive than L1):"
echo "    <2F>   =$avg_2Fl (expected: $exp_2Fl, rel. deviation: $reldev_avg_2Fl)"
echo "    <2F_H1>=$avg_2FH1l (expected: $exp_2FH1l, rel. deviation: $reldev_avg_2FH1l)"
echo "    <2F_L1>=$avg_2FL1l (expected: $exp_2FL1l, rel. deviation: $reldev_avg_2FL1l)"

if [ "$fail_avg_2Fl" -o "$fail_avg_2FH1l" -o "$fail_avg_2FL1l" ]; then
    echo "==> FAILED at tolerance=$tolerance_F"
    exit 1
else
    echo "==> OK at tolerance=$tolerance_F"
fi

echo "----------------------------------------------------------------------------------------------------"
echo "synthesizeLVstats Test2b: antenna patterns for varying detector sensitivity";
echo "----------------------------------------------------------------------------------------------------"

params_draw1=$(awk '!/%%/ && /[0-9]/' $params_file_ns | head -n 1)
alpha1=$(echo $params_draw1 | awk '{print $1}')
delta1=$(echo $params_draw1 | awk '{print $2}')
A_synth=$(echo $params_draw1 $Nsteps | awk '{print $12/(2*$28)}')
B_synth=$(echo $params_draw1 $Nsteps | awk '{print $13/(2*$28)}')
C_synth=$(echo $params_draw1 $Nsteps | awk '{print $14/(2*$28)}')
D_synth=$(echo $A_synth $B_synth $C_synth | awk '{print $1*$2-$3*$3}')
sqrtSn1=1
sqrtSn2=2
Sntot=$(echo $sqrtSn1 $sqrtSn2 | awk '{print 2/(1/($1*$1)+1/($2*$2))}')
corr1=$(echo $sqrtSn1 $Sntot | awk '{print ($1*$1)/$2}')
corr2=$(echo $sqrtSn2 $Sntot | awk '{print ($1*$1)/$2}')
A_H1_synth=$(echo $params_draw1 $corr1 $Nsteps | awk '{print $20*$28/$29}')
B_H1_synth=$(echo $params_draw1 $corr1 $Nsteps | awk '{print $21*$28/$29}')
C_H1_synth=$(echo $params_draw1 $corr1 $Nsteps | awk '{print $22*$28/$29}')
D_H1_synth=$(echo $params_draw1 $corr1 $Nsteps | awk '{print $23*$28*$28/($29*$29)}')
A_L1_synth=$(echo $params_draw1 $corr2 $Nsteps | awk '{print $24*$28/$29}')
B_L1_synth=$(echo $params_draw1 $corr2 $Nsteps | awk '{print $25*$28/$29}')
C_L1_synth=$(echo $params_draw1 $corr2 $Nsteps | awk '{print $26*$28/$29}')
D_L1_synth=$(echo $params_draw1 $corr2 $Nsteps | awk '{print $27*$28*$28/($29*$29)}')

## ----- run ComputeAntennaPattern
cap_cmdline="${cap_code} --IFOs=H1,L1 --timeStampsFile=$timestampsfile --outABCD=$outCAP2 --Alpha=$alpha1 --Delta=$delta1 --noiseSqrtShX=1,2"
echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi

A_cap=$(awk -v col=3 "$awk_print_wo_headers" $outCAP2)
B_cap=$(awk -v col=4 "$awk_print_wo_headers" $outCAP2)
C_cap=$(awk -v col=5 "$awk_print_wo_headers" $outCAP2)
D_cap=$(awk -v col=6 "$awk_print_wo_headers" $outCAP2)
A_H1_cap=$(awk -v col=7 "$awk_print_wo_headers" $outCAP2)
B_H1_cap=$(awk -v col=8 "$awk_print_wo_headers" $outCAP2)
C_H1_cap=$(awk -v col=9 "$awk_print_wo_headers" $outCAP2)
D_H1_cap=$(awk -v col=10 "$awk_print_wo_headers" $outCAP2)
A_L1_cap=$(awk -v col=11 "$awk_print_wo_headers" $outCAP2)
B_L1_cap=$(awk -v col=12 "$awk_print_wo_headers" $outCAP2)
C_L1_cap=$(awk -v col=13 "$awk_print_wo_headers" $outCAP2)
D_L1_cap=$(awk -v col=14 "$awk_print_wo_headers" $outCAP2)
reldev_A=$(echo $A_cap $A_synth | awk "$awk_reldev")
reldev_B=$(echo $B_cap $B_synth | awk "$awk_reldev")
reldev_C=$(echo $C_cap $C_synth | awk "$awk_reldev")
reldev_D=$(echo $D_cap $D_synth | awk "$awk_reldev")
reldev_A_H1=$(echo $A_H1_cap $A_H1_synth | awk "$awk_reldev")
reldev_B_H1=$(echo $B_H1_cap $B_H1_synth | awk "$awk_reldev")
reldev_C_H1=$(echo $C_H1_cap $C_H1_synth | awk "$awk_reldev")
reldev_D_H1=$(echo $D_H1_cap $D_H1_synth | awk "$awk_reldev")
reldev_A_L1=$(echo $A_L1_cap $A_L1_synth | awk "$awk_reldev")
reldev_B_L1=$(echo $B_L1_cap $B_L1_synth | awk "$awk_reldev")
reldev_C_L1=$(echo $C_L1_cap $C_L1_synth | awk "$awk_reldev")
reldev_D_L1=$(echo $D_L1_cap $D_L1_synth | awk "$awk_reldev")
fail_A=$(echo $reldev_A $tolerance_p | awk "$awk_isgtr")
fail_B=$(echo $reldev_B $tolerance_p | awk "$awk_isgtr")
fail_C=$(echo $reldev_C $tolerance_p | awk "$awk_isgtr")
fail_D=$(echo $reldev_D $tolerance_p | awk "$awk_isgtr")
fail_A_H1=$(echo $reldev_A_H1 $tolerance_p | awk "$awk_isgtr")
fail_B_H1=$(echo $reldev_B_H1 $tolerance_p | awk "$awk_isgtr")
fail_C_H1=$(echo $reldev_C_H1 $tolerance_p | awk "$awk_isgtr")
fail_D_H1=$(echo $reldev_D_H1 $tolerance_p | awk "$awk_isgtr")
fail_A_L1=$(echo $reldev_A_L1 $tolerance_p | awk "$awk_isgtr")
fail_B_L1=$(echo $reldev_B_L1 $tolerance_p | awk "$awk_isgtr")
fail_C_L1=$(echo $reldev_C_L1 $tolerance_p | awk "$awk_isgtr")
fail_D_L1=$(echo $reldev_D_L1 $tolerance_p | awk "$awk_isgtr")

echo "==> from lalapps_ComputeAntennaPattern:"
echo "<A>=$A_cap, <B>=$B_cap, <C>=$C_cap, <D>=$D_cap, <A_H1>=$A_H1_cap, <B_H1>=$B_H1_cap, <C_H1>=$C_H1_cap, <D_H1>=$D_H1_cap, <A_L1>=$A_L1_cap, <B_L1>=$B_L1_cap, <C_L1>=$C_L1_cap, <D_L1>=$D_L1_cap"
echo "    from lalapps_synthesizeLVStats, corrected with SnH1/Sn=$corr1 and SnL1/Sn=$corr2:"
echo "<A>=$A_synth, <B>=$B_synth, <C>=$C_synth, <D>=$D_synth, <A_H1>=$A_H1_synth, <B_H1>=$B_H1_synth, <C_H1>=$C_H1_synth, <D_H1>=$D_H1_synth, <A_L1>=$A_L1_synth, <B_L1>=$B_L1_synth, <C_L1>=$C_L1_synth, <D_L1>=$D_L1_synth"

if [ "$fail_A" -o "$fail_B" -o "$fail_C" -o "$fail_D" -o "$fail_A_H1" -o "$fail_B_H1" -o "$fail_C_H1" -o "$fail_D_H1" -o "$fail_A_L1" -o "$fail_B_L1" -o "$fail_C_L1" -o "$fail_D_L1" ]; then
    echo "==> FAILED at tolerance=$tolerance_p"
    exit 1
else
    echo "==> OK at tolerance=$tolerance_p"
fi


echo "----------------------------------------------------------------------------------------------------"
echo "synthesizeLVstats Test2c: LR-stats for varying detector sensitivity";
echo "----------------------------------------------------------------------------------------------------"

stats_draw1=$(awk '!/%%/ && /[0-9]/' $stats_file_ns | head -n 1)
draw1_2F=$(echo $stats_draw1 | awk '{print $5}')
draw1_2FH1=$(echo $stats_draw1 | awk '{print $6}')
draw1_2FL1=$(echo $stats_draw1 | awk '{print $7}')
draw1_LR=$(echo $stats_draw1 | awk '{print $8}')

draw1_LR_recomp=$(echo $draw1_2F $draw1_2FH1 $draw1_2FL1 $log10e $Fstar0 $pL $rH1 $rL1 | awk '{print ( 0.5*$1 - log( (1-$6)*exp($5) + $6*( exp(0.5*$2)*$7 + exp(0.5*$3)*$8 ) ) ) * $4}')
reldev_LR_draw1=$(echo $draw1_LR $draw1_LR_recomp | awk "$awk_reldev")
fail_LR_draw1=$(echo $reldev_LR_draw1 $tolerance_LR | awk "$awk_isgtr")

echo "for draw1, 2F=$draw1_2F, 2FH1=$draw1_2FH1, 2FL1=$draw1_2FL1:"
echo "           LR    =$draw1_LR  -- recomputed externally: $draw1_LR_recomp -- rel. deviation: $reldev_LR_draw1)"

if [ "$fail_LR_draw1" ]; then
    echo "==> FAILED at tolerance=$tolerance_LR"
    exit 1
else
    echo "==> OK at tolerance=$tolerance_LR"
    echo
    echo "========== OK. All synthesizeLVstats tests PASSED. =========="
    echo
fi

## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm -rf $testDir
    echo "Cleaned up."
fi

## restore original locale, just in case someone source'd this file
export LC_ALL=$LC_ALL_old

exit $retstatus
