## functions

function assert {
    if eval test "$@"; then
        echo "OK: 'test $@' passed"
    else
        echo "ERROR: 'test $@' failed"
        exit 1
    fi
}

## generate noise-only SFTs

lalapps_Makefakedata_v5 \
    --randSeed=1234 \
    --startTime=858459411 --duration=1036800 --Tsft=1800 \
    --fmin=138.0 --Band=5.0 \
    --IFOs=H1,L1 --sqrtSX=3e-23,3e-23 \
    --outSingleSFT --outSFTdir=.

## run lalapps_ComputeFstatMCUpperLimit

h0_ref=7.47583e-25
FDR_ref=0.05

lalapps_ComputeFstatMCUpperLimit \
    --alpha=6.12 --delta=1.02 \
    --freq=140.0 --freq-band=1.0 \
    --sft-patt='*.sft' \
    --loudest-2F=65 --false-dism="${FDR_ref}" --initial-h0="${h0_ref}" \
    --output-file=ul.dat

## check that MC run finished properly

line=`sed -n '$p' ul.dat`
assert "X${line}" = "X%DONE"

echo "OK: MC run finished properly"

## check that MC run converged

line=`grep -v '^%' ul.dat | sed -n '$p'`
eval ${line}
assert "X${MC_trials}" != X
assert "X${h0}" != X
assert "X${FDR_MC_int}" != X
assert "X${FDR_2F_dist}" != X

awk_script='{
  if ( ($1 - $2)^2 > (0.05 * $2)^2 ) {
    print "ERROR: relative error between " $1 " and " $2 " is too large";
    exit 1;
  } else {
    print "OK: relative error between " $1 " and " $2 " is acceptable";
  }
}'

echo "${h0} ${h0_ref}" | awk "${awk_script}"
echo "${FDR_MC_int} ${FDR_ref}" | awk "${awk_script}"
echo "${FDR_2F_dist} ${FDR_ref}" | awk "${awk_script}"

echo "OK: MC run converged"
