## create good and bad SFTs
SFTwrite

## capture output of lalapps_dumpSFT run on all good SFTs with -H/-d/-t flags,
## capture standard output (without comments which can contain Git hashes),
## check that standard error is empty
for sft in SFT-good SFT-test*; do
    for flag in H d t; do
        echo "lalapps_dumpSFT -$flag -i ./$sft"
        lalapps_dumpSFT -$flag -i ./$sft 2>stderr.txt | grep -v '^%' >stdout-$flag-$sft.txt
        if [ -s stderr.txt ]; then
            echo "ERROR: lalapps_dumpSFT -$flag -i ./$sft should not write to standard error"
            exit 1
        fi
    done
done

## create reference tarball to make updating it a bit easier
mkdir -p newtarball/
for file in stdout-*.txt; do
    cp $file newtarball/ref-$file
done
cd newtarball/
tar zcf ../new_testdumpSFT.tar.gz *
cd ..
rm -rf newtarball/

## compare standard output to reference results
for file in stdout-*.txt; do
    if ! diff -s $file ref-$file; then
        echo "ERROR: $file and ref-$file should be equal"
        exit 1
    fi
done

## run again with each flag to test --Nmax option
for flag in H d t; do
    basecl1="lalapps_dumpSFT -$flag --Nmax 1 -i ./SFT-good"
    out1="stdout-$flag-SFT-good-first.txt"
    cl1="$basecl1 | grep -v '^%' > $out1"
    echo $cl1
    eval $cl1
    if [ -s stderr.txt ]; then
        echo "ERROR: '$basecl1' should not write to standard error"
        exit 1
    fi
    Nlines=$(wc -l < $out1)
    basecl2="lalapps_dumpSFT -$flag -i ./SFT-good"
    out2="stdout-$flag-SFT-good-all.txt"
    cl2="$basecl2 | grep -v '^%' | head -n $Nlines > $out2"
    echo $cl2
    eval $cl2
    if [ -s stderr.txt ]; then
        echo "ERROR: '$basecl2' should not write to standard error"
        exit 1
    fi
    if ! diff -s $out1 $out2; then
        echo "ERROR: $out1 and $out2 should be equal"
        exit 1
    fi
done
