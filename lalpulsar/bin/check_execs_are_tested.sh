set -e

if test "${HAVE_PYTHON}" = false; then

    echo "WARNING: cannot check that all C executables under lalpulsar/bin/ have been tested"
    echo "         without Python, since some tests for C executables are written in Python"
    echo

    exit 77

fi

tested_logfiles=`find ${LAL_TEST_BUILDDIR} -name '.tested_*.log' | sort`

not_tested=`awk '
    BEGIN {
        anyprog = 0
    }
    $1 == "programs" {
        anyprog = 1
        for (i = 3; i <= NF; ++i) {
            prog[$i] = 1
            progpath[$i] = $2 "/" $i
        }
    }
    $1 == "tested" {
        for (i = 2; i <= NF; ++i) {
            tested[$i] = 1
        }
    }
    END {
        if (anyprog) {
            for (e in prog) {
                if (!tested[e]) {
                    print progpath[e]
                }
            }
        } else {
            print "NONE"
        }
    }
' ${tested_logfiles} /dev/null | sort`

if test "X${not_tested}" = XNONE; then

    echo "ERROR: test suite under lalpulsar/bin/ could not be verified"

    exit 1

fi

if test "X${not_tested}" != X; then

    echo "ERROR: the following C executables under lalpulsar/bin/ have not been tested:"
    echo
    for not_tested_exec in ${not_tested}; do
        echo ${not_tested_exec}
    done
    echo

    exit 1

fi

echo "OK: all C executables under lalpulsar/bin/ have been tested"
