# Driver script for Weave tests

# Exit as soon as any error occurs
set -e

# Skip test if requested
if test "x$1" = xskip; then
    echo "$0: skipping test '$2'"
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
srcdir=$(cd ${LAL_TEST_SRCDIR} && pwd)
builddir=$(cd ${LAL_TEST_BUILDDIR} && pwd)

# Build directories containing required tools
injdir=$(cd ${builddir}/../Injections && pwd)
sftdir=$(cd ${builddir}/../SFTTools && pwd)
fitsdir=$(cd ${builddir}/../FITSTools && pwd)
fstatdir=$(cd ${builddir}/../Fstatistic && pwd)

# Create directory for test
testdir="${builddir}/${scriptname}.testdir"
if [ -d "${testdir}" ]; then
    echo "$0: removing contents of directory '${testdir}'"
    rm -rf "${testdir}/*"
else
    mkdir -p "${testdir}"
fi

# Run test in test directory
cd "${testdir}"
echo "$0: running test '${script}' in directory '${testdir}'"
echo
source "${script}"
cd "${builddir}"

# Remove test directory, unless NOCLEANUP is set
if [ "X${NOCLEANUP}" = X ]; then
    echo "$0: removing directory '${testdir}'"
    rm -rf "${testdir}"
fi
