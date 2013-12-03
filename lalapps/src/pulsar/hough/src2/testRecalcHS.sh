#!/bin/sh

## test-script to test lalapps_RecalcHSCandidates by running it against lalapps_HierarchicalSearch (HS),
## namely by feeding it the HS output candidate-file as input and verifying that it correctly reproduces these
## candidates

## run all LALApps programs with memory debugging
export LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},memdbg"

## take user-arguments:
extra_args="$@"

builddir="./";
injectdir="../../Injections/"

## ----- user-controlled level of debug-output detail
if [ -n "$DEBUG" ]; then
    debug=${DEBUG}
else
    debug=0	## default=quiet
fi

##---------- names of codes and input/output files
code_MFD="${injectdir}lalapps_Makefakedata_v4"
code_HS="${builddir}/lalapps_HierarchicalSearch"
code_RC="${builddir}/lalapps_RecalcHSCandidates"

if [ -n "${LALPULSAR_DATADIR}" ]; then
    EEPHEM="${LALPULSAR_DATADIR}/earth00-19-DE405.dat"
    SEPHEM="${LALPULSAR_DATADIR}/sun00-19-DE405.dat"
    code_MFD="${code_MFD} -E ${LALPULSAR_DATADIR}"
    code_HS="${code_HS} --ephemE=${EEPHEM} --ephemS=${SEPHEM}"
    code_RC="${code_RC} --ephemE=${EEPHEM} --ephemS=${SEPHEM}"
fi

## ----- parameters
SFTsH1="./allSFTs_H1.sft"
SFTsL1="./allSFTs_L1.sft"

HS_OUT="./HS_out.dat"
RC_OUT1="./RC_out1.dat"
RC_OUT2="./RC_out2.dat"

Freq0=100.5
dFreq0=0.5e-2
FreqBand=1e-2
tStack=3600
nStacks=10
Tspan=36000	## should be >= nStacks * tStack
f1dot=0
f1dotBand=0
df1dot=1

echo "----- STEP 1: produce some fake data:"
cmdline="$code_MFD -I H1 --outSingleSFT --outSFTbname=$SFTsH1 -G 820108814 --duration=$Tspan --noiseSqrtSh=3e-23 --fmin=100 --Band=1 --randSeed=1"

echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cmdline' ..."
    exit 1
fi

cmdline="$code_MFD -I L1 --outSingleSFT --outSFTbname=$SFTsL1 -G 820108814 --duration=$Tspan --noiseSqrtSh=3e-23 --fmin=100 --Band=1 --randSeed=2"

echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cmdline' ..."
    exit 1
fi


patch=0.8
pixel=2

echo "----- STEP 2: '$HS_OUT' -- run original HierarchicalSearch code on this data, producing some candidates:"
cmdline="$code_HS --method=0 --DataFiles=\"$SFTsH1;$SFTsL1\" --skyRegion='(1,1)' --Freq=$Freq0 --FreqBand=$FreqBand --dFreq=$dFreq0 --f1dot=$f1dot --f1dotBand=$f1dotBand --df1dot=$df1dot --tStack=$tStack --nStacksMax=$nStacks --printCand1 --nf1dotRes=1 --semiCohPatchX=$patch --semiCohPatchY=$patch --pixelFactor=$pixel --semiCohToplist -o $HS_OUT"

echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cmdline' ..."
    exit 1
fi

echo "## ----- STEP 3a: '$RC_OUT1' -- run 'RecalcHSCandidates' on this output candidate file"
cmdline="$code_RC --DataFiles1=\"$SFTsH1;$SFTsL1\" --tStack=$tStack --nStacksMax=$nStacks --followupList=$HS_OUT --WU_Freq=$Freq0 --WU_dFreq=$dFreq0 --semiCohPatchX=$patch --semiCohPatchY=$patch --pixelFactor=$pixel --fnameout=$RC_OUT1 --outputFX=true"

echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cmdline' ..."
    exit 1
fi

echo "## ----- STEP 3b: '$RC_OUT2' -- re-run RecalcHSCandidates on its own output to test stability"
## NOTE: this time we don't use any HS-frequency bug corrections (--WU_Freq), as RC outputs frequencies correctly!
cmdline="$code_RC --DataFiles1=\"$SFTsH1;$SFTsL1\" --tStack=$tStack --nStacksMax=$nStacks --followupList=$RC_OUT1 --WU_dFreq=$dFreq0 --semiCohPatchX=$patch --semiCohPatchY=$patch --pixelFactor=$pixel --fnameout=$RC_OUT2 --outputFX=true"

echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cmdline' ..."
    exit 1
fi

## currently not really testing these results against HS, just output them for user to look at:
echo "$HS_OUT: --------------------"
sort $HS_OUT
echo "$RC_OUT1: --------------------"
sort $RC_OUT1

sort $RC_OUT1 > s1.dat
sort $RC_OUT2 > s2.dat
if ! diff -q s1.dat s2.dat; then
    echo "$RC_OUT1 and RC_OUT2 differ ... something is wrong!"
    exit 1;
else
    echo "$RC_OUT1 and RC_OUT2 are identical ... OK!"
fi

if [ -z "$NOCLEANUP" ]; then
    rm $HS_OUT $RC_OUT1 $RC_OUT2 $SFTsH1 $SFTsL1 s1.dat s2.dat
fi
