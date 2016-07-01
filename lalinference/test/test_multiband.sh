#!/usr/bin/env bash

# Exit with failure as soon as a test fails
#set -e

echo "Testing BNS: TaylorF2"
./LALInferenceMultiBandTest --psdlength 1000 --psdstart 1 --seglen 64 --srate 4096 --trigtime 0 --ifo H1 --H1-channel LALSimAdLIGO --H1-cache LALSimAdLIGO --dataseed 1324 --fix-chirpmass 1.218 --fix-q 1.0 --disable-spin --margphi --approximant TaylorF2 --0noise

echo "-------------------------------------------"
echo "Testing BNS: IMRPhenomP"
./LALInferenceMultiBandTest --psdlength 1000 --psdstart 1 --seglen 64 --srate 4096 --trigtime 0 --ifo H1 --H1-channel LALSimAdLIGO --H1-cache LALSimAdLIGO --dataseed 1324 --fix-chirpmass 1.218 --fix-q 1.0 --disable-spin --margphi --approximant IMRPhenomPv2 --0noise

echo "-------------------------------------------"
echo "Testing BBH: IMRPhenomP"
./LALInferenceMultiBandTest --psdlength 1000 --psdstart 1 --seglen 64 --srate 4096 --trigtime 0 --ifo H1 --H1-channel LALSimAdLIGO --H1-cache LALSimAdLIGO --dataseed 1324 --fix-chirpmass 10.0 --fix-q 0.7 --margphi --approximant IMRPhenomPv2 --0noise

