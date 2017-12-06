#!/bin/sh
# Shell script that runs independent test of F_+ and F_x using IndependentDetResponseTest

## allow 'make test' to work from builddir != srcdir
if [ -z "${srcdir}" ]; then
    srcdir=.
fi

echo "./IndependentDetResponseTest -c ${srcdir}/indTestDetResLHO.cfg"
./IndependentDetResponseTest -c ${srcdir}/indTestDetResLHO.cfg || exit 1

echo "./IndependentDetResponseTest -c ${srcdir}/indTestDetResLLO.cfg"
./IndependentDetResponseTest -c ${srcdir}/indTestDetResLLO.cfg || exit 1

echo "./IndependentDetResponseTest -c ${srcdir}/indTestDetResGEO.cfg"
./IndependentDetResponseTest -c ${srcdir}/indTestDetResGEO.cfg || exit 1

## all tests passed!
exit 0
