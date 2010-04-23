#!/bin/sh
# Shell script that runs independent test of F_+ and F_x using LALIndependentTestDetResponse

## allow 'make test' to work from builddir != srcdir
if [ -z "${srcdir}" ]; then
    srcdir=.
fi

echo "./LALIndependentTestDetResponse -c ${srcdir}/indTestDetResLHO.cfg"
./LALIndependentTestDetResponse -c ${srcdir}/indTestDetResLHO.cfg || exit 1

echo "./LALIndependentTestDetResponse -c ${srcdir}/indTestDetResLLO.cfg"
./LALIndependentTestDetResponse -c ${srcdir}/indTestDetResLLO.cfg || exit 1

echo "./LALIndependentTestDetResponse -c ${srcdir}/indTestDetResGEO.cfg"
./LALIndependentTestDetResponse -c ${srcdir}/indTestDetResGEO.cfg || exit 1

## all tests passed!
exit 0

