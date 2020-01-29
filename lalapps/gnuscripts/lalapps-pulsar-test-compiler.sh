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

# Test script name and location
scriptname=$(expr "X$1" : "X.*/\([^/]*\)\.sh$")
script="${LAL_TEST_SRCDIR}/${scriptname}.sh"
[ -f "${script}" ]

# Create directory for test
testdir="${LAL_TEST_BUILDDIR}/${scriptname}.testdir"
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
    ( echo *.fits | xargs -n 1 "${LAL_TEST_BUILDDIR}/../FITSTools/lalapps_fits_header_list" | grep UNCLEAN ) && exit 1
    cd "${LAL_TEST_BUILDDIR}"
else
    echo "No reference tarball ${reftarball}"
fi
echo

# Create links to any executables in lalapps/src/pulsar/ used by test script, and set PATH
for execfile in $(find "${LAL_TEST_BUILDDIR}/.." "${LAL_TEST_SRCDIR}/.." -maxdepth 2 -type f -perm -u+x); do
    execname=$(expr "X${execfile}" : "X.*/\([^/]*\)$")
    linkfile="${testdir}/${execname}"
    if grep -q "${execname}" "${script}"; then
        [ -h "${linkfile}" ] || ln -s "${execfile}" "${linkfile}"
    fi
done
export PATH="${testdir}:${PATH}"
echo "PATH=${PATH}"
echo

# List contents of test directory
echo "--- Running test in directory ${testdir} ---"
echo
ls -l "${testdir}"
echo

# Run test in test directory
echo "--- Running test $1 ---"
echo
cd "${testdir}"
export TIMEFORMAT=$'real %R\nuser %R\nsys  %R'
time bash -c "set -e; source ${script}; echo '--- Successfully ran test ${script} ---'"
cd "${LAL_TEST_BUILDDIR}"
echo

# Remove test directory, unless NOCLEANUP is set
if [ "X${NOCLEANUP}" = X ]; then
    echo "--- Removing directory ${testdir} ---"
    echo
    rm -rf "${testdir}"
fi
