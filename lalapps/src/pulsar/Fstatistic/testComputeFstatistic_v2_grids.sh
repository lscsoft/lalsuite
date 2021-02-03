##---------- names of codes and input/output files
mfdv4_code="lalapps_Makefakedata_v4"
mfdv5_code="lalapps_Makefakedata_v5"
cfsv2_code="lalapps_ComputeFstatistic_v2"
cmp_code="lalapps_compareFstats"
LTC_code="lalapps_ComputeFstatLatticeCount"

# ---------- fixed parameter of our test-signal
Tsft=1800;
startTime=711595934
duration=144000		## 40 hours

Alpha=2.0
AlphaBand=0.1
dAlpha=0.04
Delta=-0.5
DeltaBand=0.1
dDelta=0.04

h0=1
cosi=-0.3

psi=0.6
phi0=1.5

Freq=100.12345
FreqBand=0.1
dFreq=0.04
f1dot=-1e-13
f1dotBand=1e-14
df1dot=4e-15

mfd_FreqBand=4.0;
mfd_fmin=$(echo $Freq $mfd_FreqBand | awk '{printf "%g", $1 - $2 / 2.0}');

echo "mfd_fmin = $mfd_fmin"

noiseSqrtSh=0

# reference files in fe3d08fc7c81d72e9be2c4251f16346f886db3f7
# were generated with old default number of running median bins
RngMedWindow=50

## ------------------------------------------------------------

if [ "$noiseSqrtSh" != 0 ]; then
    sqrtSh=$noiseSqrtSh
    haveNoise=true;
else
    sqrtSh=1;	## for SemiAnalyticF signal-only case
    haveNoise=false;
fi

IFO=L1

##--------------------------------------------------
## test starts here
##--------------------------------------------------

echo
echo "----------------------------------------------------------------------"
echo " STEP 1: Generate Fake SFTs"
echo "----------------------------------------------------------------------"
echo

## --- for grid types 0,1,2,3,6, generate 40 hours of contiguous SFTs

## create SFT directory
SFTdir_40h="./SFTs_40h"
rm -rf $SFTdir_40h
mkdir $SFTdir_40h

# build MFD_v4 command line
mfd_CL=" --Alpha=$Alpha --Delta=$Delta --IFO=$IFO --Tsft=$Tsft --startTime=$startTime --duration=$duration --h0=$h0 --cosi=$cosi --psi=$psi --phi0=$phi0"
mfd_CL="${mfd_CL} --fmin=$mfd_fmin --Band=$mfd_FreqBand --Freq=$Freq --outSFTbname=$SFTdir_40h --f1dot=$f1dot"
if [ "$haveNoise" = true ]; then
    mfd_CL="$mfd_CL --noiseSqrtSh=$sqrtSh";
fi

## generate SFTs
cmdline="$mfdv4_code $mfd_CL";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdv4_code' ..."
    exit 1
fi

## --- for grid types 8,9, generate 2 SFTs spaced 5 days apart

## create SFT directory
SFTdir_5d="./SFTs_5d"
rm -rf $SFTdir_5d
mkdir $SFTdir_5d

## create timestamps file
cat <<EOF >"${SFTdir_5d}/timestamps.txt"
800000000 0
800432000 0
EOF

## build MFD_v5 command line
mfd_CL="--outSingleSFT=true --outSFTdir=${SFTdir_5d} --IFOs=H1 --sqrtSX=1.0 --timestampsFiles=${SFTdir_5d}/timestamps.txt --fmin=100 --Band=1.0 --randSeed=12345"

## generate SFTs
cmdline="$mfdv5_code $mfd_CL";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error: something failed when running '$mfdv5_code' ..."
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo " STEP 2: run CFS_v2 for various grid-types"
echo "----------------------------------------------------------------------"
echo

basecmds=("" "" "" "" "" "" "" "" "" "") # length 10 to cover grid types 0-9

## common arguments for grid types 0,1,2,3,6
sky_CL="--Alpha=$Alpha --AlphaBand=$AlphaBand --dAlpha=$dAlpha --Delta=$Delta --DeltaBand=$DeltaBand --dDelta=$dDelta"
spin_CL="--Freq=$Freq --FreqBand=$FreqBand --dFreq=$dFreq --f1dot=$f1dot --f1dotBand=$f1dotBand --df1dot=$df1dot"
cfs_CL="--DataFiles='${SFTdir_40h}/*.sft' --TwoFthreshold=0 --Dterms=16 --FstatMethod=DemodOptC --RngMedWindow=$RngMedWindow $extra_args"
if [ "$haveNoise" = false ]; then
    cfs_CL="$cfs_CL --assumeSqrtSX=1"
fi

## ----- grid=0 : flat grid
echo "CFSv2 using gridType=0:"
basecmds[0]="$cfsv2_code $sky_CL $spin_CL $cfs_CL --gridType=0"
cmdline="${basecmds[0]} --outputFstat=./testCFSv2_grid0.dat";
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cmdline' ..."
    exit 1
fi

## ----- grid=1 : isotropic
echo "CFSv2 using gridType=1:"
basecmds[1]="$cfsv2_code $sky_CL $spin_CL $cfs_CL --gridType=1"
cmdline="${basecmds[1]} --outputFstat=./testCFSv2_grid1.dat";
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cmdline' ..."
    exit 1
fi

## ----- grid=2 : metric grid
echo "CFSv2 using gridType=2:"
basecmds[2]="$cfsv2_code $sky_CL $spin_CL $cfs_CL --gridType=2 --metricType=1 --metricMismatch=0.1"
cmdline="${basecmds[2]} --outputFstat=./testCFSv2_grid2.dat";
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cmdline' ..."
    exit 1
fi

## ----- grid=3 : same metric grid, specified as 'grid-file' passed on the commandline
echo "CFSv2 using gridType=3:"
basecmds[3]="$cfsv2_code $sky_CL $spin_CL $cfs_CL --gridType=3 --gridFile='{2.000000953674316 -0.4704365134239197; 2.068386077880859 -0.4704365134239197; 2.099475860595703 -0.492608368396759; 2.000000953674316 -0.4204375147819519; 2.082211971282959 -0.4204375147819519}'"
cmdline="${basecmds[3]} --outputFstat=./testCFSv2_grid3.dat";
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cmdline' ..."
    exit 1
fi

## ----- grid=6 : recompute metric grid with --gridFile option
echo "recompute gridType=2 with CFSv2 using gridType=6:"
## extract a 6-column gridFile from the previous result output-file
gridFile="testCFSv2_grid2_grid.dat";
awk_extract6='{printf "%s %s %s %s %s %s\n", $1, $2, $3, $4, $5, $6 }'
sed '/^%.*/d' ./testCFSv2_grid2.dat | awk "$awk_extract6" >${gridFile}
basecmds[6]="$cfsv2_code $cfs_CL --gridType=6 --gridFile=./${gridFile}"
cmdline="${basecmds[6]} --outputFstat=./testCFSv2_grid6.dat";
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cmdline' ..."
    exit 1
fi

## --- grid=8 : lattice tiling grid, square parameter space
echo "CFSv2 using gridType=8:"
basecmds[8]="$cfsv2_code --Alpha=6.1 --Delta=1.2 --Freq=100.4 --FreqBand=5e-4 --f1dot=-1e-10 --f1dotBand=1e-10 --DataFiles='${SFTdir_5d}/*.sft' --TwoFthreshold=0 --gridType=8 --metricMismatch=0.5 --RngMedWindow=$RngMedWindow"
cmdline="${basecmds[8]} --outputFstat=./testCFSv2_grid8.dat"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cmdline' ..."
    exit 1
fi

## --- grid=9 : lattice tiling grid, age-spindown-index parameter space
echo "CFSv2 using gridType=9:"
basecmds[9]="$cfsv2_code --Alpha=6.1 --Delta=1.2 --Freq=100.4 --FreqBand=8e-5 --spindownAge=1e11 --minBraking=2 --maxBraking=5 --DataFiles='${SFTdir_5d}/*.sft' --TwoFthreshold=0 --gridType=9 --metricMismatch=0.5 --RngMedWindow=$RngMedWindow"
cmdline="${basecmds[9]} --outputFstat=./testCFSv2_grid9.dat"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cmdline' ..."
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo " STEP 3: Compare to reference results: "
echo "----------------------------------------------------------------------"
echo

# Some of the reference files in fe3d08fc7c81d72e9be2c4251f16346f886db3f7
# were generated with old --SignalOnly flag which added +4 assuming no-noise SFTs.
# This behaviour is no longer reproduced by the newer --assumeSqrtSX option,
# but to make sure this is the only (and understood) difference from the old
# reference results, we manually subtract that extra term here.
awk_subtract4='{print $1 " " $2 " " $3 " " $4 " " $5 " " $6 " " ($7-4)}'

for n in 0 1 2 3 6 8 9; do

    ## compare results
    echo "Comparing gridType=${n}:"
    ref="./testCFSv2_grid${n}.dat.ref"
    if grep -q 'SignalOnly=TRUE' $ref; then
        cmdline="awk '/^[^%]/ $awk_subtract4' $ref > ${ref}minus4.dat"
        if ! eval $cmdline; then
            echo "OUCH... using awk to subtract 4 from reference results failed."
        fi
        cmdline="$cmp_code -1 ./testCFSv2_grid${n}.dat -2 ${ref}minus4.dat";
    else
        cmdline="$cmp_code -1 ./testCFSv2_grid${n}.dat -2 $ref";
    fi
    echo $cmdline
    if ! eval $cmdline; then
        echo "OUCH... files differ. Something might be wrong..."
        exit 2
    else
        echo "OK."
    fi

done

echo
echo "----------------------------------------------------------------------"
echo " STEP 4: Compare template counts (gridType=8|9 only):"
echo "----------------------------------------------------------------------"
echo

## --- grid=8 : lattice tiling grid, square parameter space
cmdline="$LTC_code --time-span=433800 --square=6.1,0,1.2,0,100.4,5e-4,-1e-10,1e-10 --max-mismatch=0.5 --lattice=Ans --metric=spindown"
echo $cmdline
ntemplates=`$cmdline`
if [ $? -ne 0 ]; then
    echo "Error.. something failed when running '$cmdline' ..."
    exit 1
fi
ntemplates=`echo X$ntemplates | sed 's/[^0123456789]//g'`
ntemplates_ref=`grep -v '^%' ./testCFSv2_grid8.dat | wc -l | sed 's/[^0123456789]//g'`
echo "Compare template counts (gridType=8): '$ntemplates' vs '$ntemplates_ref'"
if [ "X$ntemplates" != "X$ntemplates_ref" ]; then
    echo "OUCH... template counts differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi

## --- grid=9 : lattice tiling grid, age-spindown-index parameter space
cmdline="$LTC_code --time-span=433800 --age-braking=6.1,1.2,100.4,8e-5,1e11,2,5 --max-mismatch=0.5 --lattice=Ans --metric=spindown"
echo $cmdline
ntemplates=`$cmdline`
if [ $? -ne 0 ]; then
    echo "Error.. something failed when running '$cmdline' ..."
    exit 1
fi
ntemplates=`echo X$ntemplates | sed 's/[^0123456789]//g'`
ntemplates_ref=`grep -v '^%' ./testCFSv2_grid9.dat | wc -l | sed 's/[^0123456789]//g'`
echo "Compare template counts (gridType=9): '$ntemplates' vs '$ntemplates_ref'"
if [ "X$ntemplates" != "X$ntemplates_ref" ]; then
    echo "OUCH... template counts differ. Something might be wrong..."
    exit 2
else
    echo "OK."
fi

echo
echo "----------------------------------------------------------------------"
echo " STEP 5: Check --outputGrid option:"
echo "----------------------------------------------------------------------"
echo

for n in 0 1 2 3 6 8 9; do

    # run with pure grid output (no search)
    echo "CFSv2 using gridType=$n:"
    cmdline="${basecmds[$n]} --outputGrid=./testCFSv2_grid${n}_nosearch.dat"
    echo $cmdline
    if ! eval $cmdline; then
        echo "Error.. something failed when running '$cmdline' ..."
        exit 1
    fi

    ## compare results
    echo "Comparing --outputFstat vs --outputGrid for gridType=${n}:"
    grid="./testCFSv2_grid${n}_nosearch.dat"
    gridcut="./testCFSv2_grid${n}_nosearch_cut.dat"
    sed_cmd="sed '/%/d' ${grid} > ${gridcut}"
    if ! eval $sed_cmd; then
        echo "Error.. something failed when running '$sed_cmd' ..."
        exit 1
    fi
    res="./testCFSv2_grid${n}.dat"
    rescut="./testCFSv2_grid${n}_cut.dat"
    sed_cmd="sed '/%/d' ${res} | cut -d ' ' -f 1-6 > ${rescut}"
    if ! eval $sed_cmd; then
        echo "Error.. something failed when running '$sed_cmd' ..."
        exit 1
    fi
    diff_cmd="diff -s ${gridcut} ${rescut}"
    echo $diff_cmd
    if ! eval $diff_cmd; then
        echo "error: columns 1-6 of $grid and $res should be equal."
        exit 1
    else
        echo "OK, columns 1-6 match."
    fi

done
