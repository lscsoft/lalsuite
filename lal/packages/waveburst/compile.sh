#!/usr/bin/env bash
make
pushd $LALPREFIX/src/lal/lib
pwd
make clean
cd ..
pwd
make
popd
pwd

