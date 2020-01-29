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
if [ "X${LAL_TEST_SRCDIR}" = X ]; then
  echo "LAL_TEST_SRCDIR is not set"
  exit 1
fi
if [ "X${LAL_TEST_BUILDDIR}" = X ]; then
  echo "LAL_TEST_BUILDDIR is not set"
  exit 1
fi

# Test script name and location
scriptname=$(expr "X$1" : "X.*/\([^/]*\)\.sh$")
script="${LAL_TEST_SRCDIR}/${scriptname}.sh"
if [ ! -f "${script}" ]; then
  echo "Test script '${script}' does not exist"
  exit 1
fi
if [ -x "${script}" ]; then
  echo "Test script '${script}' should not be executable"
  exit 1
fi

# Create directory for test
testdir="${LAL_TEST_BUILDDIR}/${scriptname}.testdir"
rm -rf "${testdir}"
if [ -d "${testdir}" ]; then
  echo "Could not remove test directory '${testdir}'"
  exit 1
fi
mkdir -p "${testdir}"

# Extract any reference results, and check validity
reftarball="${LAL_TEST_SRCDIR}/${scriptname}.tar.gz"
if [ -f ${reftarball} ]; then
    echo "Extracting reference tarball ${reftarball}"
    cd "${testdir}"
    tar xf ${reftarball}
    for file in $(find . -type f); do
        (
            case "${file}" in
                *.txt)
                    grep UNCLEAN "${file}"
                    ;;
                *.fits)
                    "${LAL_TEST_BUILDDIR}/../FITSTools/lalapps_fits_header_list" "${file}" | grep UNCLEAN
                    ;;
                *)
                    exit 1
                    ;;
            esac
        ) && {
            echo "File '${file}' appears to have been built from a git checkout of LALSuite with uncommitted changes, and therefore is not reproducible"
            exit 1
        }
    done
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
export TIMEFORMAT=$'\n\n\nreal %R\nuser %R\nsys  %R'
time bash -c "\
set -e; \
export LAL_TEST_SRCDIR=/dev/null/; \
export LAL_TEST_BUILDDIR=/dev/null/; \
export srcdir=/dev/null/; \
export builddir=/dev/null/; \
source '${script}'" && status=$? || status=$?
cd "${LAL_TEST_BUILDDIR}"
echo
case $status in
    0)
        echo "--- Test $1 ran successfully ---"
        ;;
    77)
        echo "--- Test $1 was skipped ---"
        ;;
    *)
        echo "--- Test $1 exited with status ${status} ---"
        NOCLEANUP=1
        ;;
esac
echo

# Remove test directory, unless NOCLEANUP is set
if [ "X${NOCLEANUP}" = X ]; then
    echo "--- Removing directory ${testdir} ---"
    echo
    rm -rf "${testdir}"
fi

# Return test script exit status
exit ${status}
