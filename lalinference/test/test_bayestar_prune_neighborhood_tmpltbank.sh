#!/bin/sh -e
# Nothing fancy: for now, just run the script and check that it does not crash.
../python/bayestar_prune_neighborhood_tmpltbank test_bayestar_sim_to_tmpltbank.xml --mass1 1.4 --mass2 1.4 --snr 8 \
    -o test_bayestar_prune_neighborhood_tmpltbank.xml
