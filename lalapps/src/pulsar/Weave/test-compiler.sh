# Driver script for Weave tests
set -e
echo "--- Test compiler is $0 ---"
echo

# Skip test if requested
if test "x$1" = xskip; then
    echo "--- Skipping test $2 ---"
    echo
    exit 77;
fi

# Source and build directories
[ "X${LAL_TEST_SRCDIR}" != X ]
[ "X${LAL_TEST_BUILDDIR}" != X ]
export srcdir=$(cd ${LAL_TEST_SRCDIR} && pwd)
export builddir=$(cd ${LAL_TEST_BUILDDIR} && pwd)

# Build directories containing required tools
export injdir=$(cd ${builddir}/../Injections && pwd)
export sftdir=$(cd ${builddir}/../SFTTools && pwd)
export fitsdir=$(cd ${builddir}/../FITSTools && pwd)

# Test script name and location
scriptname=$(expr "X$1" : "X.*/\([^/]*\)\.sh$")
script="${srcdir}/${scriptname}.sh"
[ -f "${script}" ]

# Create directory for test
testdir="${builddir}/${scriptname}.testdir"
echo "--- Running test in directory ${testdir} ---"
echo
rm -rf "${testdir}"
[ ! -d "${testdir}" ]
mkdir -p "${testdir}"

# Extract any reference results, and check validity
reftarball="${srcdir}/${scriptname}.tar.gz"
if [ -f ${reftarball} ]; then
    tar xf ${reftarball}
    cd "${testdir}"
    ( echo *.txt | xargs -n 1 cat | grep UNCLEAN ) && exit 1
    ( echo *.fits | xargs -n 1 ${fitsdir}/lalapps_fits_header_list | grep UNCLEAN ) && exit 1
    cd "${builddir}"
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
