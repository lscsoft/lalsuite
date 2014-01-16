#!/bin/sh

## set LAL debug level
echo "Setting LAL_DEBUG_LEVEL=${LAL_DEBUG_LEVEL:-msglvl1,memdbg}"
export LAL_DEBUG_LEVEL

## take user-arguments
extra_args="$@"

## allow 'make test' to work from builddir != srcdir
if [ -z "${srcdir}" ]; then
    srcdir=`dirname $0`
fi

builddir="./";
injectdir="../Injections/"
fdsdir="../FDS_isolated/"

##---------- names of codes
cap_code="${builddir}lalapps_ComputeAntennaPattern"
pds_code="${builddir}lalapps_PrintDetectorState"
mfd_code="${injectdir}lalapps_Makefakedata_v4"
pfs_code="${fdsdir}lalapps_PredictFStat"

tolerance=1e-3
Tsft=1800

## awk commands needed for testing
awk_reldev='{printf "%.2e", sqrt(($1-$2)*($1-$2))/(0.5*sqrt(($1+$2)*($1+$2))) }'
awk_isgtr='{if($1>$2) {print "1"}}'

# ---------- common test parameters
IFO=H1
timestamp1=852443819
## CAP has included offset of Tsft/2, so always need to pass that manually to PDS
timestamp1_pds=$( echo $timestamp1 $Tsft | awk '{printf "%d", $1+$2/2}' )
timestamp2=852445619
timestamp3=852447419
alpha=4.649850989853494
delta=-0.506281802989210

# ---------- temporary files
outCAP=antenna_pattern_test.dat
outPDS=detector_state_test.dat
skygridfile=./skygrid_test.dat
timestampsfile=./timestamps_test.dat
sftfile=./H1_test.sft
outPFS=./pfs_test.dat


## if a previous run failed, have to delete some files to avoid appending
## (the others are forcefully recreated by their respective lalapps)
if [ -f $outPDS ]; then
    rm $outPDS
fi
if [ -f $skygridfile ]; then
    rm $skygridfile
fi
if [ -f $timestampsfile ]; then
    rm $timestampsfile
fi

echo "----------------------------------------------------------------------------------------------------"
echo "ComputeAntennaPattern Test0: internal consistency of different computeAM implementations in PrintDetectorState";
echo "----------------------------------------------------------------------------------------------------"

pds_cmdline="${pds_code} -I $IFO -a $alpha -d $delta -t $timestamp1_pds >> $outPDS"
echo $pds_cmdline;
if ! eval $pds_cmdline; then
    echo "Error.. something failed when running '$pds_code' ..."
    exit 1
fi

pds_out_lalcomp=$(grep LALComputeAM ${outPDS} | tr ',' ' ')
a_lalcomp=$( echo $pds_out_lalcomp | awk '{print $3}' )
b_lalcomp=$( echo $pds_out_lalcomp | awk '{print $4}' )
echo "==> LALComputeAM:                    a=$a_lalcomp, b=$b_lalcomp"
pds_out_lalget=$(grep LALGetAMCoeffs ${outPDS} | tr ',' ' ')
a_lalget=$( echo $pds_out_lalget | awk '{print $3}' )
b_lalget=$( echo $pds_out_lalget | awk '{print $4}' )
echo "    LALGetAMCoeffs:                  a=$a_lalget, b=$b_lalget"
pds_out_lalnew=$(grep LALNewGetAMCoeffs ${outPDS} | tr ',' ' ')
a_lalnew=$( echo $pds_out_lalnew | awk '{print $3}' )
b_lalnew=$( echo $pds_out_lalnew | awk '{print $4}' )
echo "    LALNewGetAMCoeffs:               a=$a_lalnew, b=$b_lalnew"
pds_out_xlal=$(grep XLALComputeAntennaPatternCoeffs ${outPDS} | tr ',' ' ')
a_xlal=$( echo $pds_out_xlal | awk '{print $3}' )
b_xlal=$( echo $pds_out_xlal | awk '{print $4}' )
echo "    XLALComputeAntennaPatternCoeffs: a=$a_xlal, b=$b_xlal"

reldev_a_lalget=$(echo $a_lalcomp $a_lalget | awk "$awk_reldev")
reldev_b_lalget=$(echo $b_lalcomp $b_lalget | awk "$awk_reldev")
reldev_a_lalnew=$(echo $a_lalcomp $a_lalnew | awk "$awk_reldev")
reldev_b_lalnew=$(echo $b_lalcomp $b_lalnew | awk "$awk_reldev")
reldev_a_xlal=$(echo $a_lalcomp $a_xlal | awk "$awk_reldev")
reldev_b_xlal=$(echo $b_lalcomp $b_xlal | awk "$awk_reldev")

fail_a_get=$(echo $reldev_a_lalget $tolerance | awk "$awk_isgtr")
fail_b_get=$(echo $reldev_b_lalget $tolerance | awk "$awk_isgtr")
fail_a_new=$(echo $reldev_a_lalnew $tolerance | awk "$awk_isgtr")
fail_b_new=$(echo $reldev_b_lalnew $tolerance | awk "$awk_isgtr")
fail_a_xlal=$(echo $reldev_a_xlal $tolerance | awk "$awk_isgtr")
fail_b_xlal=$(echo $reldev_b_xlal $tolerance | awk "$awk_isgtr")

if [ "$fail_a_get" -o "$fail_b_get" -o "$fail_a_new" -o "$fail_b_new" -o "$fail_a_xlal" -o "$fail_b_xlal" ]; then
    echo "==> FAILED at tolerance=$tolerance"
    exit 1
else
    echo "==> OK at tolerance=$tolerance"
fi


echo "----------------------------------------------------------------------------------------------------"
echo "ComputeAntennaPattern Test1: single-sky-point, single-timestamp comparison to PrintDetectorState";
echo "----------------------------------------------------------------------------------------------------"

## ----- run ComputeAntennaPattern
cap_cmdline="${cap_code} --IFOs=$IFO --timeGPS=$timestamp1 --outputFile=$outCAP --Alpha=$alpha --Delta=$delta"
echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi

outCAPstripped=$(sed -e"/^%%.*/d" $outCAP)
a_cap=$(echo $outCAPstripped | awk '{print $4}')
b_cap=$(echo $outCAPstripped | awk '{print $5}')
echo "==> lalapps_ComputeAntennaPattern:   a=$a_cap, b=$b_cap"

reldev_a_cap_xlal=$(echo $a_xlal $a_cap | awk "$awk_reldev")
reldev_b_cap_xlal=$(echo $b_xlal $b_cap | awk "$awk_reldev")
fail_a_cap=$(echo $reldev_a_cap_xlal $tolerance | awk "$awk_isgtr")
fail_b_cap=$(echo $reldev_b_cap_xlal $tolerance | awk "$awk_isgtr")

if [ "$fail_a_cap" -o "$fail_b_cap" ]; then
    echo "==> FAILED at tolerance=$tolerance"
    exit 1
else
    echo "==> OK at tolerance=$tolerance"
fi


echo "----------------------------------------------------------------------------------------------------"
echo "ComputeAntennaPattern Test2: sky grid";
echo "----------------------------------------------------------------------------------------------------"

## ----- produce skygrid
alpha1=0.0
delta1=0.0
alpha2=0.0
delta2=0.5
alpha3=3.0
delta3=-0.5
printf "%s %s\n%s %s\n%s %s" "$alpha1" "$delta1" "$alpha2" "$delta2" "$alpha3" "$delta3" >> $skygridfile

## ----- run PrintDetectorState 3 times
rm $outPDS
pds_cmdline="${pds_code} -I $IFO -a $alpha1 -d $delta1 -t $timestamp1_pds >> $outPDS"
echo $pds_cmdline;
if ! eval $pds_cmdline; then
    echo "Error.. something failed when running '$pds_code' ..."
    exit 1
fi
pds_out_xlal=$(grep XLALComputeAntennaPatternCoeffs ${outPDS} | tr ',' ' ')
a1_pds=$( echo $pds_out_xlal | awk '{print $3}')
b1_pds=$( echo $pds_out_xlal | awk '{print $4}')

rm $outPDS
pds_cmdline="${pds_code} -I $IFO -a $alpha2 -d $delta2 -t $timestamp1_pds >> $outPDS"
echo $pds_cmdline;
if ! eval $pds_cmdline; then
    echo "Error.. something failed when running '$pds_code' ..."
    exit 1
fi
pds_out_xlal=$(grep XLALComputeAntennaPatternCoeffs ${outPDS} | tr ',' ' ')
a2_pds=$( echo $pds_out_xlal | awk '{print $3}')
b2_pds=$( echo $pds_out_xlal | awk '{print $4}')

rm $outPDS
pds_cmdline="${pds_code} -I $IFO -a $alpha3 -d $delta3 -t $timestamp1_pds >> $outPDS"
echo $pds_cmdline;
if ! eval $pds_cmdline; then
    echo "Error.. something failed when running '$pds_code' ..."
    exit 1
fi
pds_out_xlal=$(grep XLALComputeAntennaPatternCoeffs ${outPDS} | tr ',' ' ')
a3_pds=$( echo $pds_out_xlal | awk '{print $3}')
b3_pds=$( echo $pds_out_xlal | awk '{print $4}')

## ----- run ComputeAntennaPattern
cap_cmdline="${cap_code} --IFOs=$IFO --timeGPS=$timestamp1 --outputFile=$outCAP --skyGridFile=$skygridfile"

echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi

outCAP1=$(sed -e"/^%%.*/d" $outCAP | head -n 1) ## first line
outCAP2=$(sed -e"/^%%.*/d" $outCAP | head -n 2 | tail -n 1) ## second line
outCAP3=$(sed -e"/^%%.*/d" $outCAP | tail -n 2) ## third = next-to-last line (last one is blank)

a1_cap=$(echo $outCAP1 |  awk '{print $4}')
b1_cap=$(echo $outCAP1 |  awk '{print $5}')
a2_cap=$(echo $outCAP2 |  awk '{print $4}')
b2_cap=$(echo $outCAP2 |  awk '{print $5}')
a3_cap=$(echo $outCAP3 |  awk '{print $4}')
b3_cap=$(echo $outCAP3 |  awk '{print $5}')
echo "==> alpha=0.0 delta= 0.0: lalapps_ComputeAntennaPattern: a=$a1_cap, b=$b1_cap / lalapps_PrintDetectorState: a=$a1_pds, b=$b1_pds"
echo "    alpha=0.0 delta= 0.5: lalapps_ComputeAntennaPattern: a=$a2_cap, b=$b2_cap / lalapps_PrintDetectorState: a=$a2_pds, b=$b2_pds"
echo "    alpha=3.0 delta=-0.5: lalapps_ComputeAntennaPattern: a=$a3_cap, b=$b3_cap / lalapps_PrintDetectorState: a=$a3_pds, b=$b3_pds"

reldev_a1=$(echo $a1_pds $a1_cap | awk "$awk_reldev")
reldev_b1=$(echo $b1_pds $b1_cap | awk "$awk_reldev")
fail_a1=$(echo $reldev_a1 $tolerance | awk "$awk_isgtr")
fail_b1=$(echo $reldev_b1 $tolerance | awk "$awk_isgtr")
reldev_a2=$(echo $a2_pds $a2_cap | awk "$awk_reldev")
reldev_b2=$(echo $b2_pds $b2_cap | awk "$awk_reldev")
fail_a2=$(echo $reldev_a2 $tolerance | awk "$awk_isgtr")
fail_b2=$(echo $reldev_b2 $tolerance | awk "$awk_isgtr")
reldev_a3=$(echo $a3_pds $a3_cap | awk "$awk_reldev")
reldev_b3=$(echo $b3_pds $b3_cap | awk "$awk_reldev")
fail_a3=$(echo $reldev_a3 $tolerance | awk "$awk_isgtr")
fail_b3=$(echo $reldev_b3 $tolerance | awk "$awk_isgtr")

if [ "$fail_a1" -o "$fail_b1" -o "$fail_a2" -o "$fail_b2" -o "$fail_a3" -o "$fail_b3" ]; then
    echo "==> FAILED at tolerance=$tolerance"
    exit 1
else
    echo "==> OK at tolerance=$tolerance"
fi


echo "----------------------------------------------------------------------------------------------------"
echo "ComputeAntennaPattern Test3: matrix-element averaging over timestamps from file";
echo "----------------------------------------------------------------------------------------------------"

## ----- run ComputeAntennaPattern with single-timestamp input, output
cap_cmdline="${cap_code} --IFOs=$IFO --timeGPS=$timestamp1,$timestamp2,$timestamp3 --outputFile=$outCAP --Alpha=$alpha --Delta=$delta"
echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi

outCAP1=$(sed -e"/^%%.*/d" $outCAP | head -n 1) ## first line
outCAP2=$(sed -e"/^%%.*/d" $outCAP | head -n 2 | tail -n 1) ## second line
outCAP3=$(sed -e"/^%%.*/d" $outCAP | tail -n 2) ## third = next-to-last line (last one is blank)

A1_cap=$(echo $outCAP1 |  awk '{print $6}')
B1_cap=$(echo $outCAP1 |  awk '{print $7}')
C1_cap=$(echo $outCAP1 |  awk '{print $8}')
A2_cap=$(echo $outCAP2 |  awk '{print $6}')
B2_cap=$(echo $outCAP2 |  awk '{print $7}')
C2_cap=$(echo $outCAP2 |  awk '{print $8}')
A3_cap=$(echo $outCAP3 |  awk '{print $6}')
B3_cap=$(echo $outCAP3 |  awk '{print $7}')
C3_cap=$(echo $outCAP3 |  awk '{print $8}')

## ----- externally compute mean for test
Amean=$( echo $A1_cap $A2_cap $A3_cap | awk '{print ($1+$2+$3)/3}' )
Bmean=$( echo $B1_cap $B2_cap $B3_cap | awk '{print ($1+$2+$3)/3}' )
Cmean=$( echo $C1_cap $C2_cap $C3_cap | awk '{print ($1+$2+$3)/3}' )
Dmean=$( echo $Amean $Bmean $Cmean | awk '{print $1*$2-$3*$3}' )

## ----- make timestampsfile
printf "%s 0\n%s 0\n%s 0" "$timestamp1" "$timestamp2" "$timestamp3" >> $timestampsfile

## ----- run ComputeAntennaPattern with timestampsfile input, averaged output
cap_cmdline="${cap_code} --IFOs=$IFO --timeStampsFile=$timestampsfile --outputFile=$outCAP --Alpha=$alpha --Delta=$delta --averageABCD"
echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi

outCAPstripped=$(sed -e"/^%%.*/d" $outCAP)
A_cap_mean=$(echo $outCAPstripped | awk '{print $3}')
B_cap_mean=$(echo $outCAPstripped | awk '{print $4}')
C_cap_mean=$(echo $outCAPstripped | awk '{print $5}')
D_cap_mean=$(echo $outCAPstripped | awk '{print $6}')
reldev_Amean=$(echo $A_cap_mean $Amean | awk "$awk_reldev")
reldev_Bmean=$(echo $B_cap_mean $Bmean | awk "$awk_reldev")
reldev_Cmean=$(echo $C_cap_mean $Cmean | awk "$awk_reldev")
reldev_Dmean=$(echo $D_cap_mean $Dmean | awk "$awk_reldev")
fail_Amean=$(echo $reldev_Amean $tolerance | awk "$awk_isgtr")
fail_Bmean=$(echo $reldev_Bmean $tolerance | awk "$awk_isgtr")
fail_Cmean=$(echo $reldev_Cmean $tolerance | awk "$awk_isgtr")
fail_Dmean=$(echo $reldev_Dmean $tolerance | awk "$awk_isgtr")
echo "==> mean externally computed: <A>=$Amean, <B>=$Bmean, <C>=$Cmean, <D>=$Dmean"
echo "    mean from --averageABCD:  <A>=$A_cap_mean, <B>=$B_cap_mean, <C>=$C_cap_mean, <D>=$D_cap_mean"

if [ "$fail_Amean" -o "$fail_Bmean" -o "$fail_Cmean" -o "$fail_Dmean" ]; then
    echo "==> FAILED at tolerance=$tolerance"
    exit 1
else
    echo "==> OK at tolerance=$tolerance"
fi

echo "----------------------------------------------------------------------------------------------------"
echo "ComputeAntennaPattern Test4: comparing A,B,C,D with PredictFStat";
echo "----------------------------------------------------------------------------------------------------"

cap_cmdline="${cap_code} --IFOs=$IFO --timeStampsFile=$timestampsfile --outputFile=$outCAP --Alpha=$alpha --Delta=$delta --averageABCD"
echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi

outCAPstripped=$(sed -e"/^%%.*/d" $outCAP)
A_cap=$(echo $outCAPstripped | awk '{print $3}')
B_cap=$(echo $outCAPstripped | awk '{print $4}')
C_cap=$(echo $outCAPstripped | awk '{print $5}')
D_cap=$(echo $outCAPstripped | awk '{print $6}')

mfd_cmdline="${mfd_code} --IFO=$IFO --outSingleSFT --outSFTbname=$sftfile --fmin=59.95 --Band=0.1 --timestampsFile=$timestampsfile"
echo $mfd_cmdline;
if ! eval $mfd_cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

pfs_cmdline="${pfs_code} --IFO=$IFO --h0=1e-24 --cosi=0 --psi=0 --phi0=0 --Freq=60 --Alpha=$alpha --Delta=$delta --DataFiles=$sftfile --outputFstat=$outPFS --SignalOnly --printFstat=0"
echo $pfs_cmdline;
if ! eval $pfs_cmdline; then
    echo "Error.. something failed when running '$pfs_code' ..."
    exit 1
fi

A_pfs=$(grep 'A =' ${outPFS} | tr -d 'A =;')
B_pfs=$(grep 'B =' ${outPFS} | tr -d 'B =;')
C_pfs=$(grep 'C =' ${outPFS} | tr -d 'C =;')
D_pfs=$(grep 'D =' ${outPFS} | tr -d 'D =;')

echo "==> lalapps_ComputeAntennaPattern:   A=$A_cap, B=$B_cap, C=$C_cap, D=$D_cap"
echo "    lalapps_PredictFStat:            A=$A_pfs,   B=$B_pfs,   C=$C_pfs,   D=$D_pfs"

reldev_A=$(echo $A_cap $A_pfs | awk "$awk_reldev")
fail_A=$(echo $reldev_A $tolerance | awk "$awk_isgtr")
reldev_B=$(echo $B_cap $B_pfs | awk "$awk_reldev")
fail_B=$(echo $reldev_B $tolerance | awk "$awk_isgtr")
reldev_C=$(echo $C_cap $C_pfs | awk "$awk_reldev")
fail_C=$(echo $reldev_C $tolerance | awk "$awk_isgtr")
reldev_D=$(echo $D_cap $D_pfs | awk "$awk_reldev")
fail_D=$(echo $reldev_D $tolerance | awk "$awk_isgtr")

if [ "$fail_A" -o "$fail_B" -o "$fail_C" -o "$fail_D" ]; then
    echo "==> FAILED at tolerance=$tolerance"
    retstatus=1
else
    echo "==> OK at tolerance=$tolerance"
fi

echo "----------------------------------------------------------------------------------------------------"
echo "ComputeAntennaPattern Test5: multi-IFO";
echo "----------------------------------------------------------------------------------------------------"

A_H1_single=$A_cap
B_H1_single=$B_cap
C_H1_single=$C_cap
D_H1_single=$D_cap

cap_cmdline="${cap_code} --IFOs=L1 --timeStampsFile=$timestampsfile --outputFile=$outCAP --Alpha=$alpha --Delta=$delta --averageABCD"
echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi

outCAPstripped=$(sed -e"/^%%.*/d" $outCAP)
A_L1_single=$(echo $outCAPstripped |awk '{print $3}')
B_L1_single=$(echo $outCAPstripped |awk '{print $4}')
C_L1_single=$(echo $outCAPstripped |awk '{print $5}')
D_L1_single=$(echo $outCAPstripped |awk '{print $6}')

cap_cmdline="${cap_code} --IFOs=H1,L1 --timeStampsFile=$timestampsfile --outputFile=$outCAP --Alpha=$alpha --Delta=$delta --averageABCD"
echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi

outCAPstripped=$(sed -e"/^%%.*/d" $outCAP)
A_H1L1=$(echo $outCAPstripped |awk '{print $3}')
B_H1L1=$(echo $outCAPstripped |awk '{print $4}')
C_H1L1=$(echo $outCAPstripped |awk '{print $5}')
D_H1L1=$(echo $outCAPstripped |awk '{print $6}')
A_H1_from_multi=$(echo $outCAPstripped |awk '{print $7}')
B_H1_from_multi=$(echo $outCAPstripped |awk '{print $8}')
C_H1_from_multi=$(echo $outCAPstripped |awk '{print $9}')
D_H1_from_multi=$(echo $outCAPstripped |awk '{print $10}')
A_L1_from_multi=$(echo $outCAPstripped |awk '{print $11}')
B_L1_from_multi=$(echo $outCAPstripped |awk '{print $12}')
C_L1_from_multi=$(echo $outCAPstripped |awk '{print $13}')
D_L1_from_multi=$(echo $outCAPstripped |awk '{print $14}')

reldev_A_H1=$(echo $A_H1_single $A_H1_from_multi | awk "$awk_reldev")
fail_A_H1=$(echo $reldev_A_H1 $tolerance | awk "$awk_isgtr")
reldev_B_H1=$(echo $B_H1_single $B_H1_from_multi | awk "$awk_reldev")
fail_B_H1=$(echo $reldev_B_H1 $tolerance | awk "$awk_isgtr")
reldev_C_H1=$(echo $C_H1_single $C_H1_from_multi | awk "$awk_reldev")
fail_C_H1=$(echo $reldev_C_H1 $tolerance | awk "$awk_isgtr")
reldev_D_H1=$(echo $D_H1_single $D_H1_from_multi | awk "$awk_reldev")
fail_D_H1=$(echo $reldev_D_H1 $tolerance | awk "$awk_isgtr")
reldev_A_L1=$(echo $A_L1_single $A_L1_from_multi | awk "$awk_reldev")
fail_A_L1=$(echo $reldev_A_L1 $tolerance | awk "$awk_isgtr")
reldev_B_L1=$(echo $B_L1_single $B_L1_from_multi | awk "$awk_reldev")
fail_B_L1=$(echo $reldev_B_L1 $tolerance | awk "$awk_isgtr")
reldev_C_L1=$(echo $C_L1_single $C_L1_from_multi | awk "$awk_reldev")
fail_C_L1=$(echo $reldev_C_L1 $tolerance | awk "$awk_isgtr")
reldev_D_L1=$(echo $D_L1_single $D_L1_from_multi | awk "$awk_reldev")
fail_D_L1=$(echo $reldev_D_L1 $tolerance | awk "$awk_isgtr")

A_H1L1_sum=$(echo $A_H1_single $A_L1_single | awk '{print ($1+$2)/2}')
B_H1L1_sum=$(echo $B_H1_single $B_L1_single | awk '{print ($1+$2)/2}')
C_H1L1_sum=$(echo $C_H1_single $C_L1_single | awk '{print ($1+$2)/2}')
D_H1L1_sum=$(echo $A_H1L1_sum $B_H1L1_sum $C_H1L1_sum | awk '{print $1*$2-$3*$3}')

reldev_A_H1L1=$(echo $A_H1L1 $A_H1L1_sum | awk "$awk_reldev")
fail_A_H1L1=$(echo $reldev_A_H1L1 $tolerance | awk "$awk_isgtr")
reldev_B_H1L1=$(echo $B_H1L1 $B_H1L1_sum | awk "$awk_reldev")
fail_B_H1L1=$(echo $reldev_B_H1L1 $tolerance | awk "$awk_isgtr")
reldev_C_H1L1=$(echo $C_H1L1 $C_H1L1_sum | awk "$awk_reldev")
fail_C_H1L1=$(echo $reldev_C_H1L1 $tolerance | awk "$awk_isgtr")
reldev_D_H1L1=$(echo $D_H1L1 $D_H1L1_sum | awk "$awk_reldev")
fail_D_H1L1=$(echo $reldev_D_H1L1 $tolerance | awk "$awk_isgtr")

echo "==> H1 from single-IFO run:           A=$A_H1_single, B=$B_H1_single, C=$C_H1_single, D=$D_H1_single"
echo "    H1 from multi-IFO run:            A=$A_H1_from_multi, B=$B_H1_from_multi, C=$C_H1_from_multi, D=$D_H1_from_multi"
echo "    L1 from single-IFO run:           A=$A_L1_single, B=$B_L1_single, C=$C_L1_single, D=$D_L1_single"
echo "    L1 from multi-IFO run:            A=$A_L1_from_multi, B=$B_L1_from_multi, C=$C_L1_from_multi, D=$D_L1_from_multi"
echo "    H1L1 summed from single-IFO runs: A=$A_H1L1_sum, B=$B_H1L1_sum, C=$C_H1L1_sum, D=$D_H1L1_sum"
echo "    H1L1 from multi-IFO run:          A=$A_H1L1, B=$B_H1L1, C=$C_H1L1, D=$D_H1L1"

if [ "$fail_A_H1" -o "$fail_B_H1" -o "$fail_C_H1" -o "$fail_D_H1" -o "$fail_A_L1" -o "$fail_B_L1" -o "$fail_C_L1" -o "$fail_D_L1" -o "$fail_A_H1L1" -o "$fail_B_H1L1" -o "$fail_C_H1L1" -o "$fail_D_H1L1" ]; then
    echo "==> FAILED at tolerance=$tolerance"
    retstatus=1
else
    echo "==> OK at tolerance=$tolerance"
    echo
    echo "========== OK. All ComputeAntennaPattern tests PASSED. =========="
    echo
fi

## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm $outCAP
    rm $outPDS
    rm $skygridfile
    rm $timestampsfile
    rm $sftfile
    rm $outPFS
    echo "Cleaned up."
fi

exit $retstatus
