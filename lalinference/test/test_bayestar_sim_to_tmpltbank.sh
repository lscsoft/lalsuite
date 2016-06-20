#!/bin/sh -e
# Nothing fancy: for now, just run the script and check that it does not crash.
../python/bayestar_sim_to_tmpltbank $srcdir/HL-INJECTIONS_1_TEST-1000000000-10.xml \
    -o test_bayestar_sim_to_tmpltbank.xml
