## create good and bad SFTs
SFTwrite

## these SFTs have to pass
for sft in SFT-good SFT-test*; do
    echo "lalapps_SFTvalidate $sft"
    if ! lalapps_SFTvalidate $sft; then
        echo "lalapps_SFTvalidate failed SFT $sft; should have passed"
        exit 1
    fi
    echo
done
echo "lalapps_SFTvalidate SFT-good SFT-test*"
if ! lalapps_SFTvalidate SFT-good SFT-test*; then
    echo "lalapps_SFTvalidate failed SFTs SFT-good SFT-test*; should have passed"
    exit 1
fi
echo
echo "printf '%s\n' SFT-good '' SFT-test* | lalapps_SFTvalidate >valid-sfts.txt"
if ! printf '%s\n' SFT-good '' SFT-test* | lalapps_SFTvalidate >valid-sfts.txt; then
    echo "lalapps_SFTvalidate failed SFTs SFT-good SFT-test*; should have passed"
    exit 1
fi
echo

## these SFTs have to fail
for sft in SFT-bad*; do
    echo "lalapps_SFTvalidate $sft"
    if lalapps_SFTvalidate $sft; then
        echo "lalapps_SFTvalidate passed SFT $sft; should have failed"
        exit 1
    fi
    echo
done
echo "lalapps_SFTvalidate SFT-bad*"
if lalapps_SFTvalidate SFT-bad*; then
    echo "lalapps_SFTvalidate passed SFTs SFT-bad*; should have failed"
    exit 1
fi
echo
echo "printf '%s\n' SFT-bad* | lalapps_SFTvalidate"
if printf '%s\n' SFT-bad* | lalapps_SFTvalidate; then
    echo "lalapps_SFTvalidate passed SFTs SFT-bad*; should have failed"
    exit 1
fi
echo
echo "printf '%s\n' SFT-good SFT-bad* SFT-test* | lalapps_SFTvalidate >valid-sfts-2.txt"
if printf '%s\n' SFT-good SFT-bad* SFT-test* | lalapps_SFTvalidate >valid-sfts-2.txt; then
    echo "lalapps_SFTvalidate passed SFT-good SFT-bad* SFT-test*; should have failed"
    exit 1
fi
echo
echo "diff valid-sfts.txt valid-sfts-2.txt"
if ! diff valid-sfts.txt valid-sfts-2.txt; then
    echo "valid-sfts.txt and valid-sfts-2.txt should compare equal"
    exit 1
fi
echo
