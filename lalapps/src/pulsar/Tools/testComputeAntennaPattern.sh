#!/bin/sh

## set LAL debug level
echo "Setting LAL_DEBUG_LEVEL=${LAL_DEBUG_LEVEL:-msglvl1,memdbg}"
export LAL_DEBUG_LEVEL

LC_ALL_old=$LC_ALL
export LC_ALL=C

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
tolerance_pfs=1 ## more lenient because PFS has noise fluctuations from MFD
Tsft=1800

## awk commands needed for testing
awk_reldev='{printf "%.2e", sqrt(($1-$2)^2)/(0.5*sqrt(($1+$2)^2))}'
awk_isgtr='{if($1>$2) {print "1"}}'
awk_print_wo_headers='!/%%/ && /[0-9]/ {print $col}'
## awk line filtering: !/%%/ filters out header lines and /[0-9]/ filters out blank lines (require at least one number)
awk_avg='!/%%/ && /[0-9]/ { total += $col; count++ } END { print total/count }'
awk_avg_sq='!/%%/ && /[0-9]/ { total += $col*$col; count++ } END { print total/count }'
awk_avg_prod='!/%%/ && /[0-9]/ { total += $col1*$col2; count++ } END { print total/count }'

# ---------- common test parameters
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
timestampsfile_H1=./timestamps_test_H1.dat
timestampsfile_L1=./timestamps_test_L1.dat
sftfile_base=./sft_test_
sftfile_H1="${sftfile_base}H1"
sftfile_L1="${sftfile_base}L1"
outPFS=./pfs_test.dat


## if a previous run failed, have to delete some files to avoid appending
## (the others are forcefully recreated by their respective lalapps)
if [ -f $outPDS ]; then
    rm $outPDS
fi
if [ -f $skygridfile ]; then
    rm $skygridfile
fi
if [ -f $timestampsfile_H1 ]; then
    rm $timestampsfile_H1
fi
if [ -f $timestampsfile_L1 ]; then
    rm $timestampsfile_L1
fi

echo "----------------------------------------------------------------------------------------------------"
echo "ComputeAntennaPattern Test0: internal consistency of different computeAM implementations in PrintDetectorState";
echo "----------------------------------------------------------------------------------------------------"

pds_cmdline="${pds_code} -I H1 -a $alpha -d $delta -t $timestamp1_pds >> $outPDS"
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
echo "ComputeAntennaPattern Test1: single-sky-point, single-timestamp comparison of a(t), b(t) to PrintDetectorState";
echo "----------------------------------------------------------------------------------------------------"

## ----- run ComputeAntennaPattern
cap_cmdline="${cap_code} --IFOs=H1 --timeGPS=$timestamp1 --outab=$outCAP --Alpha=$alpha --Delta=$delta"
echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi

a_cap=$(awk -v col=4 "$awk_print_wo_headers" $outCAP)
b_cap=$(awk -v col=5 "$awk_print_wo_headers" $outCAP)
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
echo "ComputeAntennaPattern Test2: a(t), b(t) over sky grid";
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
pds_cmdline="${pds_code} -I H1 -a $alpha1 -d $delta1 -t $timestamp1_pds >> $outPDS"
echo $pds_cmdline;
if ! eval $pds_cmdline; then
    echo "Error.. something failed when running '$pds_code' ..."
    exit 1
fi
pds_out_xlal=$(grep XLALComputeAntennaPatternCoeffs ${outPDS} | tr ',' ' ')
a1_pds=$( echo $pds_out_xlal | awk '{print $3}')
b1_pds=$( echo $pds_out_xlal | awk '{print $4}')

rm $outPDS
pds_cmdline="${pds_code} -I H1 -a $alpha2 -d $delta2 -t $timestamp1_pds >> $outPDS"
echo $pds_cmdline;
if ! eval $pds_cmdline; then
    echo "Error.. something failed when running '$pds_code' ..."
    exit 1
fi
pds_out_xlal=$(grep XLALComputeAntennaPatternCoeffs ${outPDS} | tr ',' ' ')
a2_pds=$( echo $pds_out_xlal | awk '{print $3}')
b2_pds=$( echo $pds_out_xlal | awk '{print $4}')

rm $outPDS
pds_cmdline="${pds_code} -I H1 -a $alpha3 -d $delta3 -t $timestamp1_pds >> $outPDS"
echo $pds_cmdline;
if ! eval $pds_cmdline; then
    echo "Error.. something failed when running '$pds_code' ..."
    exit 1
fi
pds_out_xlal=$(grep XLALComputeAntennaPatternCoeffs ${outPDS} | tr ',' ' ')
a3_pds=$( echo $pds_out_xlal | awk '{print $3}')
b3_pds=$( echo $pds_out_xlal | awk '{print $4}')

## ----- run ComputeAntennaPattern
cap_cmdline="${cap_code} --IFOs=H1 --timeGPS=$timestamp1 --outab=$outCAP --skyGridFile=$skygridfile"

echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi

outCAP1=$(awk '!/%%/ && /[0-9]/' $outCAP | head -n 1) ## first line
outCAP2=$(awk '!/%%/ && /[0-9]/' $outCAP | head -n 2 | tail -n 1) ## second line
outCAP3=$(awk '!/%%/ && /[0-9]/' $outCAP | tail -n 3 | tail -n 1) ## third line

a1_cap=$(echo $outCAP1 |  awk '{print $4}')
b1_cap=$(echo $outCAP1 |  awk '{print $5}')
a2_cap=$(echo $outCAP2 |  awk '{print $4}')
b2_cap=$(echo $outCAP2 |  awk '{print $5}')
a3_cap=$(echo $outCAP3 |  awk '{print $4}')
b3_cap=$(echo $outCAP3 |  awk '{print $5}')

# read a1_cap a2_cap a3_cap <<< $(awk '!/%%/ && /[0-9]/ {print $4}' $outCAP)
# read b1_cap b2_cap b3_cap <<< $(awk '!/%%/ && /[0-9]/ {print $5}' $outCAP)

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
cap_cmdline="${cap_code} --IFOs=H1 --timeGPS=$timestamp1,$timestamp2,$timestamp3 --outab=$outCAP --Alpha=$alpha --Delta=$delta"
echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi

## ----- externally compute mean ABCD for test
Amean=$(awk -v col=4 "$awk_avg_sq" $outCAP)
Bmean=$(awk -v col=5 "$awk_avg_sq" $outCAP)
Cmean=$(awk -v col1=4 -v col2=5 "$awk_avg_prod" $outCAP)
Dmean=$( echo $Amean $Bmean $Cmean | awk '{print $1*$2-$3*$3}' )

## ----- make timestampsfile
printf "%s 0\n%s 0\n%s 0" "$timestamp1" "$timestamp2" "$timestamp3" >> $timestampsfile_H1

## ----- run ComputeAntennaPattern with timestampsfile input, direct average ABCD output
cap_cmdline="${cap_code} --IFOs=H1 --timeStampsFiles=$timestampsfile_H1 --outABCD=$outCAP --Alpha=$alpha --Delta=$delta"
echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi

A_cap_mean=$(awk -v col=3 "$awk_print_wo_headers" $outCAP)
B_cap_mean=$(awk -v col=4 "$awk_print_wo_headers" $outCAP)
C_cap_mean=$(awk -v col=5 "$awk_print_wo_headers" $outCAP)
D_cap_mean=$(awk -v col=6 "$awk_print_wo_headers" $outCAP)
reldev_Amean=$(echo $A_cap_mean $Amean | awk "$awk_reldev")
reldev_Bmean=$(echo $B_cap_mean $Bmean | awk "$awk_reldev")
reldev_Cmean=$(echo $C_cap_mean $Cmean | awk "$awk_reldev")
reldev_Dmean=$(echo $D_cap_mean $Dmean | awk "$awk_reldev")
fail_Amean=$(echo $reldev_Amean $tolerance | awk "$awk_isgtr")
fail_Bmean=$(echo $reldev_Bmean $tolerance | awk "$awk_isgtr")
fail_Cmean=$(echo $reldev_Cmean $tolerance | awk "$awk_isgtr")
fail_Dmean=$(echo $reldev_Dmean $tolerance | awk "$awk_isgtr")
echo "==> computed from a(t), b(t): <A>=$Amean,  <B>=$Bmean,  <C>=$Cmean,  <D>=$Dmean"
echo "    direct ABCD output:       <A>=$A_cap_mean, <B>=$B_cap_mean, <C>=$C_cap_mean, <D>=$D_cap_mean"

if [ "$fail_Amean" -o "$fail_Bmean" -o "$fail_Cmean" -o "$fail_Dmean" ]; then
    echo "==> FAILED at tolerance=$tolerance"
    exit 1
else
    echo "==> OK at tolerance=$tolerance"
fi

echo "----------------------------------------------------------------------------------------------------"
echo "ComputeAntennaPattern Test4: comparing A,B,C,D with PredictFStat";
echo "----------------------------------------------------------------------------------------------------"

cap_cmdline="${cap_code} --IFOs=H1 --timeStampsFiles=$timestampsfile_H1 --outABCD=$outCAP --Alpha=$alpha --Delta=$delta"
echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi

A_cap=$(awk -v col=3 "$awk_print_wo_headers" $outCAP)
B_cap=$(awk -v col=4 "$awk_print_wo_headers" $outCAP)
C_cap=$(awk -v col=5 "$awk_print_wo_headers" $outCAP)
D_cap=$(awk -v col=6 "$awk_print_wo_headers" $outCAP)

mfd_cmdline="${mfd_code} --IFO=H1 --outSingleSFT --outSFTbname=$sftfile_H1 --fmin=59.95 --Band=0.1 --timestampsFile=$timestampsfile_H1"
echo $mfd_cmdline;
if ! eval $mfd_cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

pfs_cmdline="${pfs_code} --IFO=H1 --h0=1e-24 --cosi=0 --psi=0 --phi0=0 --Freq=60 --Alpha=$alpha --Delta=$delta --DataFiles=$sftfile_H1 --outputFstat=$outPFS --SignalOnly --printFstat=0"
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
    exit 1
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

cap_cmdline="${cap_code} --IFOs=L1 --timeStampsFiles=$timestampsfile_H1 --outABCD=$outCAP --Alpha=$alpha --Delta=$delta"
echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi

A_L1_single=$(awk -v col=3 "$awk_print_wo_headers" $outCAP)
B_L1_single=$(awk -v col=4 "$awk_print_wo_headers" $outCAP)
C_L1_single=$(awk -v col=5 "$awk_print_wo_headers" $outCAP)
D_L1_single=$(awk -v col=6 "$awk_print_wo_headers" $outCAP)

cap_cmdline="${cap_code} --IFOs=H1,L1 --timeStampsFiles=$timestampsfile_H1 --outABCD=$outCAP --Alpha=$alpha --Delta=$delta"
echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi

A_H1L1=$(awk -v col=3 "$awk_print_wo_headers" $outCAP)
B_H1L1=$(awk -v col=4 "$awk_print_wo_headers" $outCAP)
C_H1L1=$(awk -v col=5 "$awk_print_wo_headers" $outCAP)
D_H1L1=$(awk -v col=6 "$awk_print_wo_headers" $outCAP)
A_H1_from_multi=$(awk -v col=7 "$awk_print_wo_headers" $outCAP)
B_H1_from_multi=$(awk -v col=8 "$awk_print_wo_headers" $outCAP)
C_H1_from_multi=$(awk -v col=9 "$awk_print_wo_headers" $outCAP)
D_H1_from_multi=$(awk -v col=10 "$awk_print_wo_headers" $outCAP)
A_L1_from_multi=$(awk -v col=11 "$awk_print_wo_headers" $outCAP)
B_L1_from_multi=$(awk -v col=12 "$awk_print_wo_headers" $outCAP)
C_L1_from_multi=$(awk -v col=13 "$awk_print_wo_headers" $outCAP)
D_L1_from_multi=$(awk -v col=14 "$awk_print_wo_headers" $outCAP)

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

A_H1L1_comb=$(echo $A_H1_single $A_L1_single | awk '{print ($1+$2)/2}')
B_H1L1_comb=$(echo $B_H1_single $B_L1_single | awk '{print ($1+$2)/2}')
C_H1L1_comb=$(echo $C_H1_single $C_L1_single | awk '{print ($1+$2)/2}')
D_H1L1_comb=$(echo $A_H1L1_comb $B_H1L1_comb $C_H1L1_comb | awk '{print $1*$2-$3*$3}')

reldev_A_H1L1=$(echo $A_H1L1 $A_H1L1_comb | awk "$awk_reldev")
fail_A_H1L1=$(echo $reldev_A_H1L1 $tolerance | awk "$awk_isgtr")
reldev_B_H1L1=$(echo $B_H1L1 $B_H1L1_comb | awk "$awk_reldev")
fail_B_H1L1=$(echo $reldev_B_H1L1 $tolerance | awk "$awk_isgtr")
reldev_C_H1L1=$(echo $C_H1L1 $C_H1L1_comb | awk "$awk_reldev")
fail_C_H1L1=$(echo $reldev_C_H1L1 $tolerance | awk "$awk_isgtr")
reldev_D_H1L1=$(echo $D_H1L1 $D_H1L1_comb | awk "$awk_reldev")
fail_D_H1L1=$(echo $reldev_D_H1L1 $tolerance | awk "$awk_isgtr")

echo "==> H1 from single-IFO run:             A=$A_H1_single, B=$B_H1_single, C=$C_H1_single, D=$D_H1_single"
echo "    H1 from multi-IFO run:              A=$A_H1_from_multi, B=$B_H1_from_multi, C=$C_H1_from_multi, D=$D_H1_from_multi"
echo "    L1 from single-IFO run:             A=$A_L1_single, B=$B_L1_single, C=$C_L1_single, D=$D_L1_single"
echo "    L1 from multi-IFO run:              A=$A_L1_from_multi, B=$B_L1_from_multi, C=$C_L1_from_multi, D=$D_L1_from_multi"
echo "    H1L1 combined from single-IFO runs: A=$A_H1L1_comb,  B=$B_H1L1_comb,  C=$C_H1L1_comb,  D=$D_H1L1_comb"
echo "    H1L1 from multi-IFO run:            A=$A_H1L1, B=$B_H1L1, C=$C_H1L1, D=$D_H1L1"

if [ "$fail_A_H1" -o "$fail_B_H1" -o "$fail_C_H1" -o "$fail_D_H1" -o "$fail_A_L1" -o "$fail_B_L1" -o "$fail_C_L1" -o "$fail_D_L1" -o "$fail_A_H1L1" -o "$fail_B_H1L1" -o "$fail_C_H1L1" -o "$fail_D_H1L1" ]; then
    echo "==> FAILED at tolerance=$tolerance"
    exit 1
else
    echo "==> OK at tolerance=$tolerance"
fi

echo "----------------------------------------------------------------------------------------------------"
echo "ComputeAntennaPattern Test6: varying detector sensitivity (compared with PFS)";
echo "----------------------------------------------------------------------------------------------------"

## compute harmonic mean of ShX
SqrtShH=3e-23
SqrtShL=6e-23
Sinv=$(echo $SqrtShH $SqrtShL | awk '{print 0.5 * (1/($1^2)+1/($2^2))}')

## need more timestamps to get decent statistics for PFS
rm $timestampsfile_H1
Nsteps=50
iTS=1
while [ $iTS -le $Nsteps ]; do
    timestamp_i=$(echo $timestamp1 $iTS $Tsft| awk '{print $1 + ($2 - 1) * $3}')
    printf "%s 0\n" "$timestamp_i" >> $timestampsfile_H1
    printf "%s 0\n" "$timestamp_i" >> $timestampsfile_L1
    iTS=$(($iTS + 1))
done
## get one extra timestamp in L1 to catch errors on unequal Nsteps
iTS=$(($iTS + 1))
timestamp_i=$(echo $timestamp1 $iTS $Tsft| awk '{print $1 + ($2 - 1) * $3}')
printf "%s 0\n" "$timestamp_i" >> $timestampsfile_L1

cap_cmdline="${cap_code} --IFOs=H1,L1 --timeStampsFiles=$timestampsfile_H1,$timestampsfile_L1 --outABCD=$outCAP --Alpha=$alpha --Delta=$delta --noiseSqrtShX=$SqrtShH,$SqrtShL"
echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi

A_H1L1=$(awk -v col=3 "$awk_print_wo_headers" $outCAP)
B_H1L1=$(awk -v col=4 "$awk_print_wo_headers" $outCAP)
C_H1L1=$(awk -v col=5 "$awk_print_wo_headers" $outCAP)
D_H1L1=$(awk -v col=6 "$awk_print_wo_headers" $outCAP)
A_H1=$(awk -v col=7 "$awk_print_wo_headers" $outCAP)
B_H1=$(awk -v col=8 "$awk_print_wo_headers" $outCAP)
C_H1=$(awk -v col=9 "$awk_print_wo_headers" $outCAP)
D_H1=$(awk -v col=10 "$awk_print_wo_headers" $outCAP)
A_L1=$(awk -v col=11 "$awk_print_wo_headers" $outCAP)
B_L1=$(awk -v col=12 "$awk_print_wo_headers" $outCAP)
C_L1=$(awk -v col=13 "$awk_print_wo_headers" $outCAP)
D_L1=$(awk -v col=14 "$awk_print_wo_headers" $outCAP)

mfd_cmdline="${mfd_code} --randSeed=1 --IFO=H1 --outSingleSFT --outSFTbname=$sftfile_H1 --fmin=59.95 --Band=0.1 --timestampsFile=$timestampsfile_H1 --noiseSqrtSh=$SqrtShH"
echo $mfd_cmdline;
if ! eval $mfd_cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

pfs_cmdline="${pfs_code} --IFO=H1 --h0=1e-24 --cosi=0 --psi=0 --phi0=0 --Freq=60 --Alpha=$alpha --Delta=$delta --DataFiles=$sftfile_H1 --outputFstat=$outPFS --printFstat=0"
echo $pfs_cmdline;
if ! eval $pfs_cmdline; then
    echo "Error.. something failed when running '$pfs_code' ..."
    exit 1
fi

A_H1_pfs=$(grep 'A =' ${outPFS} | tr -d 'A =;')
B_H1_pfs=$(grep 'B =' ${outPFS} | tr -d 'B =;')
C_H1_pfs=$(grep 'C =' ${outPFS} | tr -d 'C =;')
D_H1_pfs=$(grep 'D =' ${outPFS} | tr -d 'D =;')

reldev_A_H1=$(echo $A_H1 $A_H1_pfs | awk "$awk_reldev")
fail_A_H1=$(echo $reldev_A_H1 $tolerance_pfs | awk "$awk_isgtr")
reldev_B_H1=$(echo $B_H1 $B_H1_pfs | awk "$awk_reldev")
fail_B_H1=$(echo $reldev_B_H1 $tolerance_pfs | awk "$awk_isgtr")
reldev_C_H1=$(echo $C_H1 $C_H1_pfs | awk "$awk_reldev")
fail_C_H1=$(echo $reldev_C_H1 $tolerance_pfs | awk "$awk_isgtr")
reldev_D_H1=$(echo $D_H1 $D_H1_pfs | awk "$awk_reldev")
fail_D_H1=$(echo $reldev_D_H1 $tolerance_pfs | awk "$awk_isgtr")

mfd_cmdline="${mfd_code} --randSeed=2 --IFO=L1 --outSingleSFT --outSFTbname=$sftfile_L1 --fmin=59.95 --Band=0.1 --timestampsFile=$timestampsfile_L1 --noiseSqrtSh=$SqrtShL"
echo $mfd_cmdline;
if ! eval $mfd_cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

pfs_cmdline="${pfs_code} --IFO=L1 --h0=1e-24 --cosi=0 --psi=0 --phi0=0 --Freq=60 --Alpha=$alpha --Delta=$delta --DataFiles=$sftfile_L1 --outputFstat=$outPFS --printFstat=0"
echo $pfs_cmdline;
if ! eval $pfs_cmdline; then
    echo "Error.. something failed when running '$pfs_code' ..."
    exit 1
fi

A_L1_pfs=$(grep 'A =' ${outPFS} | tr -d 'A =;')
B_L1_pfs=$(grep 'B =' ${outPFS} | tr -d 'B =;')
C_L1_pfs=$(grep 'C =' ${outPFS} | tr -d 'C =;')
D_L1_pfs=$(grep 'D =' ${outPFS} | tr -d 'D =;')

reldev_A_L1=$(echo $A_L1 $A_L1_pfs | awk "$awk_reldev")
fail_A_L1=$(echo $reldev_A_L1 $tolerance_pfs | awk "$awk_isgtr")
reldev_B_L1=$(echo $B_L1 $B_L1_pfs | awk "$awk_reldev")
fail_B_L1=$(echo $reldev_B_L1 $tolerance_pfs | awk "$awk_isgtr")
reldev_C_L1=$(echo $C_L1 $C_L1_pfs | awk "$awk_reldev")
fail_C_L1=$(echo $reldev_C_L1 $tolerance_pfs | awk "$awk_isgtr")
reldev_D_L1=$(echo $D_L1 $D_L1_pfs | awk "$awk_reldev")
fail_D_L1=$(echo $reldev_D_L1 $tolerance_pfs | awk "$awk_isgtr")

pfs_cmdline="${pfs_code} --h0=1e-24 --cosi=0 --psi=0 --phi0=0 --Freq=60 --Alpha=$alpha --Delta=$delta --DataFiles=$sftfile_base* --outputFstat=$outPFS --printFstat=0"
echo $pfs_cmdline;
if ! eval $pfs_cmdline; then
    echo "Error.. something failed when running '$pfs_code' ..."
    exit 1
fi

A_H1L1_pfs=$(grep 'A =' ${outPFS} | tr -d 'A =;')
B_H1L1_pfs=$(grep 'B =' ${outPFS} | tr -d 'B =;')
C_H1L1_pfs=$(grep 'C =' ${outPFS} | tr -d 'C =;')
D_H1L1_pfs=$(grep 'D =' ${outPFS} | tr -d 'D =;')

reldev_A_H1L1=$(echo $A_H1L1 $A_H1L1_pfs | awk "$awk_reldev")
fail_A_H1L1=$(echo $reldev_A_H1L1 $tolerance_pfs | awk "$awk_isgtr")
reldev_B_H1L1=$(echo $B_H1L1 $B_H1L1_pfs | awk "$awk_reldev")
fail_B_H1L1=$(echo $reldev_B_H1L1 $tolerance_pfs | awk "$awk_isgtr")
reldev_C_H1L1=$(echo $C_H1L1 $C_H1L1_pfs | awk "$awk_reldev")
fail_C_H1L1=$(echo $reldev_C_H1L1 $tolerance_pfs | awk "$awk_isgtr")
reldev_D_H1L1=$(echo $D_H1L1 $D_H1L1_pfs | awk "$awk_reldev")
fail_D_H1L1=$(echo $reldev_D_H1L1 $tolerance_pfs | awk "$awk_isgtr")

echo "==> H1 CAP:   A=$A_H1, B=$B_H1, C=$C_H1, D=$D_H1"
echo "    H1 PFS:   A=$A_H1_pfs, B=$B_H1_pfs, C=$C_H1_pfs, D=$D_H1_pfs"
echo "    L1 CAP:   A=$A_L1, B=$B_L1, C=$C_L1, D=$D_L1"
echo "    L1 PFS:   A=$A_L1_pfs, B=$B_L1_pfs, C=$C_L1_pfs, D=$D_L1_pfs"
echo "    H1L1 CAP: A=$A_H1L1, B=$B_H1L1, C=$C_H1L1, D=$D_H1L1"
echo "    H1L1 PFS: A=$A_H1L1_pfs,   B=$B_H1L1_pfs,   C=$C_H1L1_pfs,   D=$D_H1L1_pfs"

if [ "$fail_A_H1" -o "$fail_B_H1" -o "$fail_C_H1" -o "$fail_D_H1" -o "$fail_A_L1" -o "$fail_B_L1" -o "$fail_C_L1" -o "$fail_D_L1" -o "$fail_A_H1L1" -o "$fail_B_H1L1" -o "$fail_C_H1L1" -o "$fail_D_H1L1" ]; then
    echo "==> FAILED at tolerance=$tolerance_pfs"
    exit 1
#     retstatus=1
else
    echo "==> OK at tolerance=$tolerance_pfs"
    echo
    echo "========== OK. All ComputeAntennaPattern tests PASSED. =========="
    echo
fi

## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm $outCAP
    rm $outPDS
    rm $skygridfile
    rm $timestampsfile_H1 $timestampsfile_L1
    rm $sftfile_H1 $sftfile_L1
    rm $outPFS
    echo "Cleaned up."
fi

## restore original locale, just in case someone source'd this file
export LC_ALL=$LC_ALL_old

exit $retstatus
