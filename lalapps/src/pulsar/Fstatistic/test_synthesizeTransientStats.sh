## functions

awk_reldev='{printf "%.2e", sqrt(($1-$2)*($1-$2))/(0.5*($1+$2)) }'
awk_isgtr='{if($1>$2) {print "1"}}'

function check_average {
    local file="$1"
    local col="$2"
    local colname="$3"
    local ref="$4"
    local awk_script='
      /^%/ { next }
      { tot += $'"$col"'; count += 1 }
      END {
        avg = tot / count;
        ref = '"$ref"';
        if ( (avg - ref)^2 > (0.05 * ref)^2 ) {
          print "ERROR: '"file=${file} col=${colname}"': relative error between values " avg " and " ref " is too large";
          exit 1;
        } else {
          print "OK: '"file=${file} col=${colname}"': relative error between values " avg " and " ref " is acceptable";
        }
      }
    '
    awk "${awk_script}" "${file}"
}

function twoF_from_atoms(){
    local file="$1"
    local awk_script='
      /^%/ { next }
      {
        A += $2;
        B += $3;
        C += $4;
        Fa_re += $5;
        Fa_im += $6;
        Fb_re += $7;
        Fb_im += $8;
        count += 1
      }
      END {
        Dinv = 1. / ( A * B - C * C );
        twoF = 2. * Dinv * ( B * ( Fa_re * Fa_re + Fa_im * Fa_im ) + A * ( Fb_re * Fb_re + Fb_im * Fb_im ) - 2.0 * C * ( Fa_re * Fb_re + Fa_im * Fb_im ) );
        print twoF
      }
    '
    twoF=$(awk "${awk_script}" "${file}")
    echo $twoF
}

## common arguments

timestamp1=818845553
duration=86400
numDraws=1000
common_args="--IFOs=H1 --dataStartGPS=$timestamp1 --dataDuration=$duration --numDraws=$numDraws --randSeed=1"

## lalapps_synthesizeTransientStats with a single square window should give the same results as lalapps_synthesizeLVStats

synthLVcode="lalapps_synthesizeLVStats"
synthTScode="lalapps_synthesizeTransientStats"

transient_args="--injectWindow-type=rect --injectWindow-tauDays=1 --injectWindow-tauDaysBand=0 --injectWindow-t0Days=0 --injectWindow-t0DaysBand=0 --searchWindow-type=rect --searchWindow-tauDays=1 --searchWindow-tauDaysBand=0 --searchWindow-t0Days=0 --searchWindow-t0DaysBand=0"

CLsynthLV="$synthLVcode ${common_args} --computeBSGL=false"
CLsynthTS="$synthTScode ${common_args} ${transient_args}"

# --- first try with Gaussian noise

cmdline="${CLsynthLV} --fixedSNR=0 --outputStats=synthLV_stats_H1_g.dat"
echo $cmdline
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$synthLVcode' ..."
    exit 1
fi
cmdline="${CLsynthTS} --fixedSNR=0 --outputStats=synthTS_stats_H1_g.dat"
echo $cmdline
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$synthTScode' ..."
    exit 1
fi

check_average "synthLV_stats_H1_g.dat" 5 "2F" 4.0
check_average "synthTS_stats_H1_g.dat" 9 "2F" 4.0

# --- now try with signals

cmdline="${CLsynthLV} --fixedSNR=4 --outputStats=synthLV_stats_H1_s.dat"
echo $cmdline
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$synthLVcode' ..."
    exit 1
fi
cmdline="${CLsynthTS} --fixedSNR=4 --outputStats=synthTS_stats_H1_s.dat"
echo $cmdline
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$synthTScode' ..."
    exit 1
fi

check_average "synthLV_stats_H1_s.dat" 5 "2F" 20.0
check_average "synthTS_stats_H1_s.dat" 9 "2F" 20.0

# --- now try multi-IFO mode

common_args="--IFOs=H1,L1 --dataStartGPS=$timestamp1 --dataDuration=$duration --numDraws=$numDraws --randSeed=1"
cmdline="$synthTScode ${common_args} ${transient_args} --fixedSNR=4 --outputStats=synthTS_stats_H1L1_s.dat"
echo $cmdline
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$synthTScode' ..."
    exit 1
fi
check_average "synthTS_stats_H1L1_s.dat" 9 "2F" 20.0

# --- multi-IFO and timestampsfiles

tsH1="ts_H1.txt"
tsL1="ts_L1.txt"
TSFT=1800
Nsfts=48
Tend=$(($timestamp1 + $duration - $TSFT)) # last SFT should *end* on first+duration
for i in `seq $timestamp1 $TSFT $Tend`; do
    echo "$i 0" >> ${tsH1}
done
for i in `seq $timestamp1 $TSFT $Tend`; do
    echo "$(($i + 900)) 0" >> ${tsL1}
done

common_args="--IFOs=H1,L1 --timestampsFiles=$tsH1,$tsL1 --numDraws=$numDraws --randSeed=1"
cmdline="$synthTScode ${common_args} ${transient_args} --fixedSNR=4 --outputStats=synthTS_stats_H1L1_s_TS.dat"
echo $cmdline
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$synthTScode' ..."
    exit 1
fi
check_average "synthTS_stats_H1L1_s_TS.dat" 9 "2F" 20.0

# test single- and multi-IFO atoms on a single draw
# very strong signal but avoiding SignalOnly with its custom +-4 handling
statsH1="synthTS_stats_H1_s_N1.dat"
statsH1L1="synthTS_stats_H1L1_s_N1.dat"
atomsH1="synthTS_atoms_H1_s"
atomsH1L1="synthTS_atoms_H1L1_s_TS"
snr=4000
common_args="--numDraws=1 --randSeed=1 --fixedSNR=$snr --computeFtotal"
cmdline="$synthTScode ${common_args} ${transient_args} --IFOs=H1 --timestampsFiles=$tsH1 --outputStats=$statsH1 --outputAtoms=$atomsH1"
echo $cmdline
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$synthTScode' ..."
    exit 1
fi
cmdline="$synthTScode ${common_args} ${transient_args} --IFOs=H1,L1 --timestampsFiles=$tsH1,$tsL1 --outputStats=$statsH1L1 --outputAtoms=$atomsH1L1"
echo $cmdline
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$synthTScode' ..."
    exit 1
fi

echo "Checking number of timestamps in atoms output..."
atomsH1="${atomsH1}_0001_of_0001.dat"
atomsH1L1="${atomsH1L1}_0001_of_0001.dat"

NatomsH1=$(grep -v % $atomsH1 | wc -l)
NatomsH1L1=$(grep -v % $atomsH1L1 | wc -l)
if [ $NatomsH1 -ne $Nsfts ]; then
    echo "Error.. number of atoms $NatomsH1 in H1 output file does not match expectation $Nsfts ..."
    exit 1
fi
Natoms2=$(echo $Nsfts | awk '{printf "%d", 2*$1}')
if [ $NatomsH1L1 -ne $Natoms2 ]; then
    echo "Error.. number of atoms $NatomsH1L1 in H1L1 output file does not match expectation 2*$Nsfts=$Natoms2 ..."
    exit 1
fi

echo "...got $NatomsH1 atoms in H1 output file and $NatomsH1L1 in H1L1 output file, as expected."

# compare the TOTAL F-stat against summing up the atoms ourselves
# (in this example, we expect maxTwoF=twoF, but in general deciding how many atoms to sum up for maxTwoF would be more work)
# and also compare recovered SNRs from sqrt(2F-4) against injected value
twoF_H1=$(grep -v % $statsH1 | tr -s " " | sed "s/^ *//g" | cut -d " " -f 6)
snr_H1=$(echo $twoF_H1 | awk '{printf "%.6f", sqrt($1-4)}')
reldev0=$(echo $snr_H1 $snr | awk "$awk_reldev")
echo "2F from stats file (H1 only): $twoF_H1 (snr=$snr_H1 off from injected $snr by $reldev0)"
twoF_H1L1=$(grep -v % $statsH1L1 | tr -s " " | sed "s/^ *//g" | cut -d " " -f 6)
snr_H1L1=$(echo $twoF_H1L1 | awk '{printf "%.6f", sqrt($1-4)}')
reldev1=$(echo $snr_H1L1 $snr | awk "$awk_reldev")
echo "2F from stats file (H1L1): $twoF_H1L1 (snr=$snr_H1L1 off from injected $snr by $reldev1)"

twoF_H1_from_atoms=$(twoF_from_atoms ${atomsH1})
reldev2=$(echo $twoF_H1_from_atoms $twoF_H1 | awk "$awk_reldev")
echo "2F from atoms (H1 only): $twoF_H1_from_atoms (off from stats file by $reldev2)"
twoF_H1L1_from_atoms=$(twoF_from_atoms ${atomsH1L1})
reldev3=$(echo $twoF_H1L1_from_atoms $twoF_H1L1 | awk "$awk_reldev")
echo "2F from atoms (H1L1): $twoF_H1L1_from_atoms (off from stats file by $reldev3)"

grep -v % $atomsH1L1 | head -n 48 > "${atomsH1L1}_firsthalf"
twoF_H1_from_H1L1_atoms=$(twoF_from_atoms "${atomsH1L1}_firsthalf")
echo "2F from atoms (H1 only from first half of H1L1 file): $twoF_H1_from_H1L1_atoms"
grep -v % $atomsH1L1 | tail -n 48 > "${atomsH1L1}_secondhalf"
twoF_L1_from_H1L1_atoms=$(twoF_from_atoms "${atomsH1L1}_secondhalf")
echo "2F from atoms (L1 only from second half of H1L1 file): $twoF_L1_from_H1L1_atoms"
snr_summed=$(echo $twoF_H1_from_H1L1_atoms $twoF_L1_from_H1L1_atoms | awk '{printf "%.6f", sqrt($1-4+$2-4) }')
reldev4=$(echo $snr_summed $snr | awk "$awk_reldev")
echo "sqrt(sum(snr2)): $snr_summed (off from injected $snr by $reldev4)"

tolerance=1e-3
echo "Checking at tolerance $tolerance..."
fail0=$(echo $reldev0 $tolerance | awk "$awk_isgtr")
fail1=$(echo $reldev1 $tolerance | awk "$awk_isgtr")
fail2=$(echo $reldev2 $tolerance | awk "$awk_isgtr")
fail3=$(echo $reldev3 $tolerance | awk "$awk_isgtr")
fail4=$(echo $reldev4 $tolerance | awk "$awk_isgtr")
if [ "$fail0" -o "$fail1" -o "$fail2" -o "$fail3" -o "$fail4" ]; then
    echo " ==> *FAILED*"
    exit 1
else
    echo " ==> OK"
fi
