## create good and bad SFTs
SFTwrite >/dev/null

goodSFTs="SFT-good SFT-test*"
badSFTs="SFT-bad*"

## these SFTs have to pass
for i in $goodSFTs; do
    echo "lalapps_SFTvalidate ./$i"
    if ! lalapps_SFTvalidate ./$i; then
        echo "lalapps_SFTvalidate failed SFT $i; should have passed"
        exit 1
    fi
    echo
done

## these SFTs have to fail
for i in $badSFTs; do
    echo "lalapps_SFTvalidate ./$i"
    if lalapps_SFTvalidate ./$i; then
        echo "lalapps_SFTvalidate passed SFT $i; should have failed"
	exit 1
    fi
    echo
done
