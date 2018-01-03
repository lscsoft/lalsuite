#!/usr/bin/env bash

# Exit with failure as soon as a test fails
#set -e

# TF2 test disabled because reference phase is inconsistent between std and frequency series version of approximant in LALSim

echo "Testing BNS: TaylorF2"
./LALInferenceMultiBandTest --psdlength 1000 --psdstart 1 --seglen 64 --srate 4096 --trigtime 0 --ifo H1 --H1-channel LALSimAdLIGO --H1-cache LALSimAdLIGO --dataseed 1324 --fix-chirpmass 1.218 --fix-q 1.0 --disable-spin --approximant TaylorF2 --0noise --amporder 0 --H1-flow 30

echo "-------------------------------------------"
echo "Testing BNS: IMRPhenomP"
./LALInferenceMultiBandTest --psdlength 1000 --psdstart 1 --seglen 64 --srate 4096 --trigtime 0 --ifo H1 --H1-channel LALSimAdLIGO --H1-cache LALSimAdLIGO --dataseed 1324 --fix-chirpmass 1.218 --fix-q 1.0 --disable-spin --approximant IMRPhenomPv2 --0noise --amporder 0 --H1-flow 30

echo "-------------------------------------------"
echo "Testing BBH: IMRPhenomP"
./LALInferenceMultiBandTest --psdlength 1000 --psdstart 1 --seglen 64 --srate 4096 --trigtime 0 --ifo H1 --H1-channel LALSimAdLIGO --H1-cache LALSimAdLIGO --dataseed 1324 --fix-chirpmass 10.0 --fix-q 0.7  --approximant IMRPhenomPv2 --0noise --amporder 0 --H1-flow 30


