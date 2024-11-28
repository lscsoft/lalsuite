set -ex

test "$(./lal_tconvert --rfc-2822 800000000)" == "Fri, 13 May 2005 06:13:07 +0000"
test "$(./lal_tconvert Fri, 13 May 2005 06:13:07 +0000)" -eq 800000000

if [[ "$(uname -m)" != "i686" ]]; then  # not 32-bit CPU
test "$(./lal_tconvert 1831518865)" == "Tue Jan 19 03:14:07 UTC 2038"   # Y2038 check
test "$(./lal_tconvert 1831518866)" == "Tue Jan 19 03:14:08 UTC 2038"   # Y2038 check
fi
