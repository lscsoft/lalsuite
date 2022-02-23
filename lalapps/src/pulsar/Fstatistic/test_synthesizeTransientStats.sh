## functions

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

## common argments

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
Tend=$(($timestamp1 + $duration))
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
