#!/bin/sh

# check for lal_simd_detect
lal_simd_detect="${LAL_TEST_BUILDDIR}/../../bin/lal_simd_detect"
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
simd_compiler=`echo ${simd_compiler}`
echo "$0: compiler supports: ${simd_compiler}"

# get instruction sets supported by machine
simd_machine=`${lal_simd_detect} | sed -n 4p`
simd_machine=`echo ${simd_machine}`
echo "$0: machine supports: ${simd_machine}"

# get common instruction sets supported by compiler and machine
simd_common=
for simd in ${simd_compiler}; do
    case " ${simd_machine} " in
        *" ${simd} "*)
            simd_common="${simd_common} ${simd}"
            ;;
        *)
            ;;
    esac
done
simd_common=`echo ${simd_common}`
echo "$0: compiler/machine supports: ${simd_common}"

# find the most advanced SSE/AVX-family instruction sets for testing
simd_test_sse=
simd_test_avx=
for simd in ${simd_common}; do
    case "${simd}" in
        SSE*)
            simd_test_sse="${simd}"
            ;;
        AVX*)
            simd_test_avx="${simd}"
            ;;
        *)
            ;;
    esac
done
simd_test="${simd_test_sse} ${simd_test_avx}"
echo "$0: testing instruction sets: ${simd_test}"

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
