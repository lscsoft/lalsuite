#!/bin/sh -e
# Nothing fancy: for now, just run the script and check that it does not crash.
../python/bayestar_littlehope $srcdir/HL-INJECTIONS_1_TEST-1000000000-10.xml \
    --template-bank test_bayestar_sim_to_tmpltbank.xml \
    --reference-psd test_bayestar_sample_model_psd.xml \
    --detector H1 --detector L1 --waveform TaylorF2threePointFivePN \
    -o test_bayestar_littlehope.xml
