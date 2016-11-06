# Driver script for Weave tests
set -e
echo "--- Test compiler is $0 ---"

# Skip test if requested
if test "x$1" = xskip; then
    echo "--- Skipping test $2 ---"
    exit 77;
fi

# Test script name and location
scriptname=$(expr "X$1" : "X.*/\([^/]*\)\.sh$")
scriptdir=$(cd $(expr "X$1" : "X\(.*\)/[^/]*$") && pwd)
script="${scriptdir}/${scriptname}.sh"
[ -f ${script} ]

# Source and build directories
[ "X${LAL_TEST_SRCDIR}" != X ]
[ "X${LAL_TEST_BUILDDIR}" != X ]
export srcdir=$(cd ${LAL_TEST_SRCDIR} && pwd)
export builddir=$(cd ${LAL_TEST_BUILDDIR} && pwd)

# Build directories containing required tools
export injdir=$(cd ${builddir}/../Injections && pwd)
export sftdir=$(cd ${builddir}/../SFTTools && pwd)
export fitsdir=$(cd ${builddir}/../FITSTools && pwd)
export fstatdir=$(cd ${builddir}/../Fstatistic && pwd)

# Create directory for test
testdir="${builddir}/${scriptname}.testdir"
if [ -d "${testdir}" ]; then
    echo "--- Removing contents of directory ${testdir} ---"
    rm -rf "${testdir}/*"
else
    mkdir -p "${testdir}"
fi

# Run test in test directory
echo "--- Running test in directory ${testdir} ---"
cd "${testdir}"
echo "--- Running test ${script} ---"
echo
time -p ${SHELL} -c "set -e; source ${script}; echo '--- Successfully ran test ${script} ---'"
cd "${builddir}"

# Remove test directory, unless NOCLEANUP is set
if [ "X${NOCLEANUP}" = X ]; then
    echo "--- Removing directory ${testdir} ---"
    rm -rf "${testdir}"
fi
