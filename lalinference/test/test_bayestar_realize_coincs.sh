#!/bin/sh -e
# Nothing fancy: for now, just run the script and check that it does not crash.
../python/bayestar_realize_coincs --detector H1 L1 \
    --reference-psd test_bayestar_sample_model_psd.xml \
    $srcdir/HL-INJECTIONS_1_TEST-1000000000-10.xml \
    -o test_bayestar_realize_coincs.xml
