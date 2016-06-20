#!/bin/sh -e
# Nothing fancy: for now, just run the script and check that it does not crash.
bayestar_sample_model_psd \
    --H1=aLIGOZeroDetHighPower --L1=aLIGOZeroDetHighPower \
    -o test_bayestar_sample_model_psd.xml
