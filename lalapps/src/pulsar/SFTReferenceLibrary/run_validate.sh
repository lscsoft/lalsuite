#!/bin/sh

## run all LALApps programs with memory debugging
export LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},memdbg"

prog=./lalapps_SFTvalidate

./lalapps_SFTwrite >/dev/null

goodSFTs="SFT-good SFT-test*"
badSFTs="SFT-bad*"

for i in $goodSFTs; do
    $prog $i
done

## this SFTs have to yield error-results
for i in $badSFTs; do
    res=`$prog $i`
    if [ ! -n $res ]; then
	exit 1;
    fi
done
