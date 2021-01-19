##---------- names of codes and input/output files
clean_code="lalapps_SFTclean"
mfd_code="lalapps_Makefakedata_v4"
dump_code="lalapps_dumpSFT"

retstatus=0

# ---------- fixed parameters of our test SFTs - injecting a very strong line
sqrtSh=1e-22
h0=5e-22
cosi=0
startTime=828002611
fmin=50.0
fBand=0.01
fmax=$(echo "$fmin + $fBand" | bc)
linefreq="50.005"
TSFT=1800
duration=$TSFT

# ---------- define comparison functions
awk_reim_abs='{printf "%.6e", sqrt($1*$1+$2*$2) }'
awk_is_much_less='{if($1<0.1*$2) {print "1"}}'

# ---------- set output file names
original_SFT="H-1_H1_${TSFT}SFT_mfdv4-${startTime}-${TSFT}.sft"
cleaned_SFT="H-1_H1_${TSFT}SFT_cleaned-${startTime}-${duration}.sft"
dump_orig="${original_SFT}.dump"
dump_clean="${cleaned_SFT}.dump"

echo "----------------------------------------------------------------------"
echo "STEP 1: creating test data with line feature..."
echo "----------------------------------------------------------------------"
echo

cmdline="${mfd_code} --IFO=H1 --outSFTbname=. --startTime=$startTime --duration=$duration --fmin=$fmin --Band=$fBand --noiseSqrtSh=$sqrtSh --lineFeature=1 --h0=$h0 --cosi=$cosi --Freq=$linefreq --randSeed=1"

echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

## ----- dump original SFT
cmdline="${dump_code} --dataOnly --SFTfiles=$original_SFT > $dump_orig"

echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$dump_code' ..."
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo "STEP 2: 'clean' the SFT without passing a linefile (should do nothing)..."
echo "----------------------------------------------------------------------"
echo

cmdline="${clean_code} --sftDir=$original_SFT --outDir=. --fMin $fmin --fMax $fmax"

echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$clean_code' ..."
    exit 1
fi

## ----- dump the "cleaned" SFT
cmdline="${dump_code} --dataOnly --SFTfiles=$cleaned_SFT > $dump_clean"

echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$dump_code' ..."
    exit 1
fi

diff_cmd="diff -s $dump_orig $dump_clean"
echo $diff_cmd
if ! eval $diff_cmd; then
    echo "error: $dump_orig and $dump_clean should be equal as long as no --linefiles are given."
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo "STEP 3: actually clean the SFT with some simple linefiles..."
echo "----------------------------------------------------------------------"
echo

# function to de-spaghettify multiple calls with different linefiles
test_cleaning () {

    linefile_content=$1
    nexpected=$2

    linefile="H1lines.txt"
    echo "writing line definition '$linefile_content' to file $linefile ..."
    echo -e $linefile_content > $linefile

    cmdline="${clean_code} --sftDir=$original_SFT --outDir=. --fMin $fmin --fMax $fmax --linefiles=$linefile"

    echo $cmdline
    if ! eval $cmdline; then
        echo "Error.. something failed when running '$clean_code' ..."
        exit 1
    fi

    ## ----- dump the cleaned SFT
    cmdline="${dump_code} --dataOnly --SFTfiles=$cleaned_SFT > $dump_clean"

    echo $cmdline
    if ! eval $cmdline; then
        echo "Error.. something failed when running '$dump_code' ..."
        exit 1
    fi

    diff_cmd="diff --suppress-common-lines -y -s $dump_orig $dump_clean"
    echo $diff_cmd
    if ! eval $diff_cmd; then
        # This diff will returned a non-0 exit status, but that's ok, we just wanted to print it. So pass.
        :
    fi
    diff_cmd="diff --suppress-common-lines -y $dump_orig $dump_clean"
    ndiff=$($diff_cmd | wc -l)
    if [ $ndiff == $nexpected ]; then
        echo "OK: $dump_orig and $dump_clean differ by exactly $nexpected line(s), as expected."
    else
        echo "error: $dump_orig and $dump_clean should differ by exactly $nexpected line(s)."
        exit 1
    fi

}

# in this first case, we actually clean where there was an injected line,
# so here we do the additional check if the absolute SFT value has really reduced
test_cleaning "50.005 0.0 1 0.0 0.0 /*single_bin_line*/" "1"
diff_tr=$($diff_cmd | tr -s ' ')
freq1=$(echo $diff_tr | cut -d ' ' -f 1)
re1=$(echo $diff_tr | cut -d ' ' -f 2)
im1=$(echo $diff_tr | cut -d ' ' -f 3)
# field 4: the "|" separator from diff -y
freq2=$(echo $diff_tr | cut -d ' ' -f 5)
re2=$(echo $diff_tr | cut -d ' ' -f 6)
im2=$(echo $diff_tr | cut -d ' ' -f 7)
if [ $freq1 != $freq2 ]; then
    echo "error: frequencies $freq1 and $freq2 differ in diff line."
    exit 1
fi
abs1=$(echo $re1 $im1 | awk "$awk_reim_abs")
abs2=$(echo $re2 $im2 | awk "$awk_reim_abs")
checked=$(echo $abs2 $abs1 | awk "$awk_is_much_less")
if [ "$checked" == "1" ]; then
    echo "OK: Strong line at $freq1 Hz has indeed been cleaned, reducing SFT data from $abs1 to $abs2."
else
    echo "error: Strong line at $freq1 Hz should have been cleaned, but SFT data is $abs2 not at least 10 times smaller than $abs1."
    exit 1
fi

# for more complicated linefile cases, we just check the number of changed bins,
# not the actual effect on the data

echo
test_cleaning "50.005 0.0 1 0.000555556 0.000555556 /*line_with_wings*/" "3"

echo
test_cleaning "50.002222222 0.002222222 4 0.0 0.0  /*comb*/" "4"

echo
test_cleaning "50.005 0.0 1 0.0 0.0 /*line1*/\n50.01 0.0 1 0.0 0.0 /*line2*/" "2"

exit $retstatus
