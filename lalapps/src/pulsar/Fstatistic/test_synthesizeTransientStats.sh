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
common_args="--IFO=H1 --dataStartGPS=$timestamp1 --dataDuration=$duration --numDraws=$numDraws --randSeed=1"

## lalapps_synthesizeTransientStats with a single square window should give the same results as lalapps_synthesizeLVStats

synthLV="lalapps_synthesizeLVStats ${common_args} --computeBSGL=false"
synthTS="lalapps_synthesizeTransientStats ${common_args} --injectWindow-type=rect --injectWindow-tauDays=1 --injectWindow-tauDaysBand=0 --injectWindow-t0Days=0 --injectWindow-t0DaysBand=0 --searchWindow-type=rect --searchWindow-tauDays=1 --searchWindow-tauDaysBand=0 --searchWindow-t0Days=0 --searchWindow-t0DaysBand=0"

# --- first try with Gaussian noise

${synthLV} --fixedSNR=0 --outputStats="synthLV_stats_g.dat"
${synthTS} --fixedSNR=0 --outputStats="synthTS_stats_g.dat"

check_average "synthLV_stats_g.dat" 5 "2F" 4.0
check_average "synthTS_stats_g.dat" 9 "2F" 4.0

# --- now try with signals

${synthLV} --fixedSNR=4 --outputStats="synthLV_stats_s.dat"
${synthTS} --fixedSNR=4 --outputStats="synthTS_stats_s.dat"

check_average "synthLV_stats_s.dat" 5 "2F" 20.0
check_average "synthTS_stats_s.dat" 9 "2F" 20.0
