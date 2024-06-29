set -ex

# doesn't need to be a real frame file, only name is parsed
touch Z-ilwd_test_frame-600000000-1.gwf

./lal_cache Z-ilwd_test_frame-600000000-1.gwf

test $(./lal_cache Z-ilwd_test_frame-600000000-1.gwf | wc -l) -eq 1

rm -f Z-ilwd_test_frame-600000000-1.gwf
