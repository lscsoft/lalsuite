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

## common arguments to lalapps_synthesizeBstatMC

common_args="--A=0.154 --B=0.234 --C=-0.0104 --E=0 --numDraws=1e4"

## generate noise samples

output_file="stats_noise.txt"
lalapps_synthesizeBstatMC ${common_args} --outputStats="${output_file}"

check_average "${output_file}"  1 "h0Nat" 0.0
check_average "${output_file}"  9 "rho2"  0.0
check_average "${output_file}" 10 "lnL"   0.0
check_average "${output_file}" 11 "2F"    4.0
check_average "${output_file}" 12 "Bstat" 3.4
echo

## generate samples at SNR=4, cosi=0, psi=0

output_file="stats_SNR-4_cosi-0.txt"
lalapps_synthesizeBstatMC ${common_args} --SNR=4 --cosi=0 --psi=0 --outputStats="${output_file}"

check_average "${output_file}"  1 "h0Nat" 1.0
check_average "${output_file}"  9 "rho2"  16.0
check_average "${output_file}" 10 "lnL"   8.0
check_average "${output_file}" 11 "2F"    20.0
check_average "${output_file}" 12 "Bstat" 8.1
echo

## generate samples at SNR=4, cosi=0, psi=0

output_file="stats_SNR-4_cosi-1.txt"
lalapps_synthesizeBstatMC ${common_args} --SNR=4 --cosi=0.99 --psi=0 --outputStats="${output_file}"

check_average "${output_file}"  1 "h0Nat" 1.0
check_average "${output_file}"  9 "rho2"  16.0
check_average "${output_file}" 10 "lnL"   8.0
check_average "${output_file}" 11 "2F"    20.0
check_average "${output_file}" 12 "Bstat" 9.9
echo

## generate samples at SNR=4, cosi random, psi-random

output_file="stats_SNR-4_cosi-random.txt"
lalapps_synthesizeBstatMC ${common_args} --SNR=4 --outputStats="${output_file}"

check_average "${output_file}"  1 "h0Nat" 1.0
check_average "${output_file}"  9 "rho2"  16.0
check_average "${output_file}" 10 "lnL"   8.0
check_average "${output_file}" 11 "2F"    20.0
check_average "${output_file}" 12 "Bstat" 9.6
echo

## generate samples at h0=10 sqrt(Sn), cosi random, psi-random

output_file="stats_h0Nat-10_cosi-random.txt"
lalapps_synthesizeBstatMC ${common_args} --h0Nat=10 --outputStats="${output_file}"

check_average "${output_file}"  1 "h0Nat" 10.0
check_average "${output_file}"  9 "rho2"  16.0
check_average "${output_file}" 10 "lnL"   8.0
check_average "${output_file}" 11 "2F"    20.0
check_average "${output_file}" 12 "Bstat" 9.8
echo
