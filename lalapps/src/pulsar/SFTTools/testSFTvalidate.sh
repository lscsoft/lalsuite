## create good and bad SFTs
SFTwrite

## these SFTs have to pass
for sft in SFT-good SFT-test*; do
    echo "lalapps_SFTvalidate ./$sft"
    if ! lalapps_SFTvalidate ./$sft; then
        echo "lalapps_SFTvalidate failed SFT $sft; should have passed"
        exit 1
    fi
    echo
done

## these SFTs have to fail
for sft in SFT-bad*; do
    echo "lalapps_SFTvalidate ./$sft"
    if lalapps_SFTvalidate ./$sft; then
        echo "lalapps_SFTvalidate passed SFT $sft; should have failed"
        exit 1
    fi
    echo
done
