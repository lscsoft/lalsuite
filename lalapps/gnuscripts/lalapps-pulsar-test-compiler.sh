# Driver script for tests of LALApps CW codes in lalapps/src/pulsar/
set -e
echo "--- Test compiler is $0 ---"
echo

# Parse command line
script_extn="$1"
echo "script_extn=${script_extn}"
shift
skip_tests="$1"
echo "skip_tests=${skip_tests}"
shift
flags=""
while [ "X$2" != X ]; do
    flags="${flags} '$1'"
    shift
done
test="$1"
echo "flags=${flags}"
echo "test=${test}"
echo

# Check for required environment variables
if [ "X${LAL_TEST_SRCDIR}" = X ]; then
  echo "LAL_TEST_SRCDIR is not set"
  exit 1
fi
if [ "X${LAL_TEST_BUILDDIR}" = X ]; then
  echo "LAL_TEST_BUILDDIR is not set"
  exit 1
fi

# Test script name, location, and extension
script_name=$(expr "X${test}" : "X.*/\(test[^/]*\)\.${script_extn}$")
if [ "X${script_name}" = X ]; then
  echo "Could not parse name of test script from '${test}'"
  exit 1
fi
script="${LAL_TEST_SRCDIR}/${script_name}.${script_extn}"
if [ ! -f "${script}" ]; then
  echo "Test script '${script}' does not exist"
  exit 1
fi
if [ -x "${script}" ]; then
  echo "Test script '${script}' should not be executable"
  exit 1
fi
case "${script_extn}" in
    sh)
        cmdline="time bash ${flags} -c 'set -e; source ${script}'"
        ;;
    py)
        cmdline="time ${PYTHON} ${flags} '${script}'"
        ;;
    *)
        echo "Test script '${script}' does not have a recognised extension"
        exit 1
        ;;
esac

# Skip test if requested
case " ${skip_tests} " in
    *" ${script_name}.${script_extn} "*)
        echo "Skipping test ${test}"
        exit 77
        ;;
esac

# Create directory for test
testdir="${LAL_TEST_BUILDDIR}/${script_name}.testdir"
rm -rf "${testdir}"
if [ -d "${testdir}" ]; then
  echo "Could not remove test directory '${testdir}'"
  exit 1
fi
mkdir -p "${testdir}"

# Extract any reference results, and check validity
reftarball="${LAL_TEST_SRCDIR}/${script_name}.tar.gz"
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
echo "--- Running test ${test}: ${cmdline} ---"
echo
cd "${testdir}"
set +e
(
    export LAL_TEST_SRCDIR=/dev/null/
    export LAL_TEST_BUILDDIR=/dev/null/
    export srcdir=/dev/null/
    export builddir=/dev/null/
    eval ${cmdline}
)
status=$?
set -e
cd "${LAL_TEST_BUILDDIR}"
echo
case $status in
    0)
        echo "--- Test ${test} ran successfully ---"
        ;;
    77)
        echo "--- Test ${test} was skipped ---"
        ;;
    *)
        echo "--- Test ${test} exited with status ${status} ---"
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
