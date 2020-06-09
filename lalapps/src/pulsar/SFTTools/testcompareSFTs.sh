## create good and bad SFTs
SFTwrite

## SFTs should compare equal
tol=1e-15
for sft in SFT-good SFT-test*; do
    echo "lalapps_compareSFTs -V -e $tol -1 ./$sft -2 ./$sft"
    if ! lalapps_compareSFTs -V -e $tol -1 ./$sft -2 ./$sft; then
        echo "lalapps_compareSFTs -V failed to compare $sft equal to itself within tolerance of $tol"
        exit 1
    fi
    echo
done

## concatenated SFTs should compare equal
tol=1e-15
cat SFT-test[123] > SFT-goodB
echo "lalapps_compareSFTs -V -e $tol -1 SFT-good -2 SFT-goodB"
if ! lalapps_compareSFTs -V -e $tol -1 SFT-good -2 SFT-goodB; then
    echo "lalapps_compareSFTs -V failed to compare SFT-good equal to SFT-goodB within tolerance of $tol"
    exit 1
fi
echo

## SFT-test1 and SFT-test8 should be very close numerically
tol=1e-11
echo "lalapps_compareSFTs -V -e $tol -1 SFT-test1 -2 SFT-test8"
if ! lalapps_compareSFTs -V -e $tol -1 SFT-test1 -2 SFT-test8; then
    echo "lalapps_compareSFTs -V failed to compare SFT-test1 equal to SFT-test8 within tolerance of $tol"
    exit 1
fi
echo
