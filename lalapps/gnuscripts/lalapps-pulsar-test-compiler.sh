# Driver script for tests of LALApps CW codes in lalapps/src/pulsar/
set -e
echo "--- Test compiler is $0 ---"
echo

# Skip test if requested
if test "x$1" = xskip; then
    echo "--- Skipping test $2 ---"
    echo
    exit 77;
fi

# Check for required environment variables
test "X${LAL_TEST_SRCDIR}" != X
test "X${LAL_TEST_BUILDDIR}" != X

# Build directories containing required tools
export builddir=$(cd ${LAL_TEST_BUILDDIR} && pwd)
export injdir=$(cd ${builddir}/../Injections && pwd)
export sftdir=$(cd ${builddir}/../SFTTools && pwd)
export fitsdir=$(cd ${builddir}/../FITSTools && pwd)

# Test script name and location
scriptname=$(expr "X$1" : "X.*/\([^/]*\)\.sh$")
script="${LAL_TEST_SRCDIR}/${scriptname}.sh"
[ -f "${script}" ]

# Create directory for test
testdir="${builddir}/${scriptname}.testdir"
echo "--- Running test in directory ${testdir} ---"
echo
rm -rf "${testdir}"
[ ! -d "${testdir}" ]
mkdir -p "${testdir}"

# Extract any reference results, and check validity
reftarball="${LAL_TEST_SRCDIR}/${scriptname}.tar.gz"
if [ -f ${reftarball} ]; then
    echo "Extracting reference tarball ${reftarball}"
    cd "${testdir}"
    tar xf ${reftarball}
    ( echo *.txt | xargs -n 1 cat | grep UNCLEAN ) && exit 1
    ( echo *.fits | xargs -n 1 ${fitsdir}/lalapps_fits_header_list | grep UNCLEAN ) && exit 1
    cd "${builddir}"
else
    echo "No reference tarball ${reftarball}"
fi
ls -l "${testdir}"
echo

# Run test in test directory
echo "--- Running test ${script} ---"
echo
cd "${testdir}"
export TIMEFORMAT=$'real %R\nuser %R\nsys  %R'
time bash -c "set -e; source ${script}; echo '--- Successfully ran test ${script} ---'"
echo
cd "${builddir}"

# Remove test directory, unless NOCLEANUP is set
if [ "X${NOCLEANUP}" = X ]; then
    echo "--- Removing directory ${testdir} ---"
    echo
    rm -rf "${testdir}"
fi
