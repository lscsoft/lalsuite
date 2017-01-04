#!/bin/sh

# check for lal_simd_detect
lal_simd_detect="${LAL_TEST_BUILDDIR}/../../src/lal_simd_detect"
if test ! -x ${lal_simd_detect}; then
    echo "$0: could not execute ${lal_simd_detect}" >&2
    exit 1
fi

# check for VectorMathTest
vector_math_test="${LAL_TEST_BUILDDIR}/VectorMathTest"
if test ! -x ${vector_math_test}; then
    echo "$0: could not execute ${vector_math_test}" >&2
    exit 1
fi

# get instruction sets supported by compiler
simd_compiler=`${lal_simd_detect} | sed -n 2p`
echo "$0: compiler supports ${simd_compiler}"

# get instruction sets supported by machine
simd_machine=`${lal_simd_detect} | sed -n 4p`
echo "$0: machine supports ${simd_machine}"

# try to test these instruction sets
simd_test="SSE AVX"

for simd in ${simd_test}; do

    # test for compiler support
    case " ${simd_compiler} " in
        *"${simd}"*)
            ;;
        *)
            echo "$0: compiler does not support ${simd}"
            continue;;
    esac

    # test for machine support
    case " ${simd_machine} " in
        *"${simd}"*)
            ;;
        *)
            echo "$0: this machine does not support ${simd}"
            continue;;
    esac

    # run VectorMathTest with this instruction sed
    LAL_SIMD_ISET="${simd}"
    export LAL_SIMD_ISET
    echo
    echo "========== testing ${simd} instruction set =========="
    ${vector_math_test} || exit 1
    echo "---------- testing ${simd} instruction set ----------"
    echo

done
