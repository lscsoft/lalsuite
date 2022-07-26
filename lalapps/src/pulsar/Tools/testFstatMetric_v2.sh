## call lalapps_FstatMetric_v2

lalapps_FstatMetric_v2 \
    --Alpha=1.42 --Delta=0.73 \
    --Freq=101.2 --f1dot=1e-9 \
    --startTime=1100000000 --duration=259200 \
    --h0=1e-25 --cosi=0.1 --psi=0.22 --phi=3.4 \
    --metricType=2 \
    --outputMetric=out.dat

## check that the expected output fields are present

fields=`awk '/%/ { next } $2 == "=" { print $1 }' out.dat | sort | xargs echo`

fields_ref="A Alpha B C D Delta DopplerCoordinates Fisher_ab Nseg Tspan aCross aPlus cosi detectorWeights detectors fkdot gDN_ij gF_ij gFav_ij gN_ij g_ij h0 m1_ij m2_ij m3_ij maxrelerr_gF maxrelerr_gPh phi0 psi refTime rho2 segmentList startTime"

if [ "X${fields}" != "X${fields_ref}" ]; then
    echo "ERROR: missing fields in out.dat"
    echo "Expected:"
    echo "   ${fields_ref}"
    echo "Got:"
    echo "   ${fields}"
    exit 1
fi
echo "OK: fields in out.dat"

## check that reported relative errors are small
if ! awk '/%/ { next } $2 == "=" && $1 ~ "err" { if ( ($3*1) > 1e-10 ) { print "ERROR: " $1 " error is too large"; exit 1; } }' out.dat; then
    exit 1
fi
echo "OK: relative errors in out.dat"
