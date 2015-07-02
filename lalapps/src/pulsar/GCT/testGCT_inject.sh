#!/bin/bash

## set LAL debug level
echo "Setting LAL_DEBUG_LEVEL=${LAL_DEBUG_LEVEL:-msglvl1,memdbg}"
export LAL_DEBUG_LEVEL

## make sure we work in 'C' locale here to avoid awk sillyness
LC_ALL_old=$LC_ALL
export LC_ALL=C

# The only thing where 'dirsep' can and should be used is in paths of the SFT files,
# as in fact SFTfileIO is the only code that requires it to be set properly. Other
# file references should be handled by the shell (or wine) and converted if necessary.
dirsep="/"

builddir="./";
injectdir="../Injections/"

if [ "`echo $1 | sed 's%.*/%%'`" = "wine" ]; then
    builddir="./";
    injectdir="$1 ./"
    fdsdir="$1 ./"
    dirsep='\'
fi

##---------- names of codes and input/output files
mfd_code="${injectdir}lalapps_Makefakedata_v5"
if test $# -eq 0 ; then
    gct_code="${builddir}lalapps_HierarchSearchGCT"
else
    gct_code="$@"
fi

testDirBase="testGCT_inject.d"
testDir="./${testDirBase}";
if [ -d "$testDir" ]; then
    rm -rf $testDir
fi
mkdir -p "$testDir"

SFTdir="${testDirBase}"

## ---------- fixed parameter of our test-signal -------------
Alpha=3.1
Delta=-0.5
h0=1.0
cosi=-0.3
psi=0.6
phi0=1.5
Freq=100.123456789
f1dot=-1e-9
f2dot=1e-18
refTime=862999869

injectionSources="{Freq=$Freq; f1dot=$f1dot; f2dot=$f2dot; Alpha=$Alpha; Delta=$Delta; psi=$psi; phi0=$phi0; h0=$h0; cosi=$cosi; refTime=$refTime}"


## perfectly targeted search in sky
AlphaSearch=$Alpha
DeltaSearch=$Delta

## Produce skygrid file for the search
skygridfile="${testDir}/tmpskygridfile.dat"
echo "$AlphaSearch $DeltaSearch" > $skygridfile

mfd_FreqBand=0.20;
mfd_fmin=100;

Dterms=8
RngMedWindow=50

gct_FreqBand=0.01
gct_F1dotBand=2.0e-10
gct_F2dotBand=0
gct_dFreq=0.000002 #"2.0e-6"
gct_dF1dot=1.0e-10
gct_dF2dot=0
gct_nCands=10

sqrtSh=1

## --------- Generate fake data set time stamps -------------
echo "----------------------------------------------------------------------"
echo " STEP 0: Generating time-stamps and segments files "
echo "----------------------------------------------------------------------"
echo

Tsft=1800
startTime=852443819
Tsegment=90000
if [ -n "$NSEGMENTS" ]; then
    Nsegments=${NSEGMENTS}
else
    Nsegments=3
fi


seggap=$(echo ${Tsegment} | awk '{printf "%.0f", $1 * 1.12345}')

tsFile_H1="${testDir}/timestampsH1.dat"  # for makefakedata
tsFile_L1="${testDir}/timestampsL1.dat"  # for makefakedata
segFile="${testDir}/segments.dat"

tmpTime=$startTime
iSeg=1
while [ $iSeg -le $Nsegments ]; do
    t0=$tmpTime
    t1=$(($t0 + $Tsegment))
    TspanHours=`echo $Tsegment | awk '{printf "%.7f", $1 / 3600.0 }'`
        ## first and last segment will be single-IFO only
    if [ $iSeg -eq 1 -o $iSeg -eq $Nsegments ]; then
        NSFT=`echo $Tsegment $Tsft |  awk '{print int(1.0 * $1 / $2 + 0.5) }'`
    else	## while all other segments are 2-IFO
        NSFT=`echo $Tsegment $Tsft |  awk '{print int(2.0 * $1 / $2 + 0.5) }'`
    fi
    echo "$t0 $t1 $TspanHours $NSFT" >> $segFile

    segs[$iSeg]=$tmpTime # save seg's beginning for later use
    echo "Segment: $iSeg of $Nsegments	GPS start time: ${segs[$iSeg]}"

    Tseg=$Tsft
    while [ $Tseg -le $Tsegment ]; do
        if [ $iSeg -ne 1 ]; then
	    echo "${tmpTime} 0" >> $tsFile_H1
        fi
        if [ $iSeg -ne $Nsegments ]; then	            ## we skip segment N for L1
	    echo "${tmpTime} 0" >> $tsFile_L1
        fi
	tmpTime=$(($tmpTime + $Tsft))
	Tseg=$(($Tseg + $Tsft))
    done

    tmpTime=$(($tmpTime + $seggap))
    iSeg=$(($iSeg + 1))
done

# ---------- common MFD commandline ----------
mfd_CL_common=" --fmin=$mfd_fmin --Band=${mfd_FreqBand} --Tsft=$Tsft --randSeed=1000 --outSingleSFT --IFOs=H1,L1 --timestampsFiles=${tsFile_H1},${tsFile_L1} --sqrtSX=${sqrtSh},${sqrtSh}"

## ---------- common GCT commandline ----------------------------------------
gct_CL_common="--gridType1=3 --nCand1=$gct_nCands --skyRegion='allsky' --Freq=$Freq --skyGridFile='$skygridfile' --printCand1 --semiCohToplist --df1dot=$gct_dF1dot --f1dot=$f1dot --f1dotBand=$gct_F1dotBand  --df2dot=$gct_dF2dot --f2dot=$f2dot --f2dotBand=$gct_F2dotBand --dFreq=$gct_dFreq --FreqBand=$gct_FreqBand --refTime=$refTime --segmentList=$segFile --Dterms=$Dterms --blocksRngMed=$RngMedWindow --FstatMethod=ResampBest --computeBSGL"

echo
echo "----------------------------------------------------------------------------------------------------"
echo " STEP 1: generate SFTs with injected signal, run HierarchSearchGCT on them "
echo "----------------------------------------------------------------------------------------------------"
echo

## generate sfts containing noise + signal
label="noise+signal"
cmdline="$mfd_code $mfd_CL_common --outSFTdir=${SFTdir} --outLabel='$label' --injectionSources=\"$injectionSources\""
echo "$cmdline";
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

## run GCT code on those sfts
rm -f checkpoint.cpt # delete checkpoint to start correctly
outfile_GCT1="${testDir}/GCT1.dat"

cmdline="$gct_code $gct_CL_common --DataFiles=\"$SFTdir$dirsep*$label*.sft\" --fnameout='$outfile_GCT1'"
echo $cmdline
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$gct_code' ..."
    exit 1
fi

topline=$(sort -nr -k7,7 $outfile_GCT1 | head -1)
freqGCT1=$(echo $topline | awk '{print $1}')
resGCT1=$(echo $topline  | awk '{print $7}')
resGCT1_H1=$(echo $topline  | awk '{print $9}')
resGCT1_L1=$(echo $topline  | awk '{print $10}')

echo
echo "----------------------------------------------------------------------------------------------------"
echo " STEP 2: generate SFTs with noise only, run HierarchSearchGCT on them with on-the-fly signal-injection"
echo "----------------------------------------------------------------------------------------------------"
echo

## generate sfts containing noise + signal
label="noiseOnly"
cmdline="$mfd_code $mfd_CL_common --outSFTdir=${SFTdir} --outLabel='$label'"
echo "$cmdline";
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

## run GCT code on those sfts
rm -f checkpoint.cpt # delete checkpoint to start correctly
outfile_GCT2="${testDir}/GCT2.dat"

cmdline="$gct_code $gct_CL_common --DataFiles=\"$SFTdir$dirsep*$label*.sft\" --fnameout='$outfile_GCT2' --injectionSources=\"$injectionSources\""
echo $cmdline
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$gct_code' ..."
    exit 1
fi

topline=$(sort -nr -k7,7 $outfile_GCT2 | head -1)
freqGCT2=$(echo $topline | awk '{print $1}')
resGCT2=$(echo $topline  | awk '{print $7}')
resGCT2_H1=$(echo $topline  | awk '{print $9}')
resGCT2_L1=$(echo $topline  | awk '{print $10}')


## ---------- compute relative differences and check against tolerance --------------------
awk_reldev='{printf "%.2e", sqrt(($1-$2)*($1-$2))/(0.5*($1+$2)) }'

freqreldev=$(echo $freqGCT1 $freqGCT2 | awk "$awk_reldev")
reldev=$(echo $resGCT1 $resGCT2 | awk "$awk_reldev")
reldev_H1=$(echo $resGCT1_H1 $resGCT2_H1 | awk "$awk_reldev")
reldev_L1=$(echo $resGCT1_L1 $resGCT2_L1 | awk "$awk_reldev")


# ---------- Check relative deviations against tolerance, report results ----------
Tolerance=5e-2	## 5%
awk_isgtr='{if($1>$2) {print "1"}}'

fail1=$(echo $freqreldev $Tolerance | awk "$awk_isgtr")
fail2=$(echo $reldev $Tolerance     | awk "$awk_isgtr")
fail3=$(echo $reldev_H1 $Tolerance     | awk "$awk_isgtr")
fail4=$(echo $reldev_L1 $Tolerance     | awk "$awk_isgtr")
if [ "$fail1" -o "$fail2" -o "$fail3" -o "$fail4" ]; then
    echo " ==> *FAILED*"
    exit 1
else
    echo " ==> OK"
fi

echo "----------------------------------------------------------------------"

## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm -rf $testDir
    echo "Cleaned up."
fi

## restore original locale, just in case someone source'd this file
export LC_ALL=$LC_ALL_old
