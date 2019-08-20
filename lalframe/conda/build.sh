#!/bin/bash

set -e

./configure \
	--prefix="${PREFIX}" \
	--enable-swig-iface \
	--disable-swig-octave \
	--disable-swig-python \
	--disable-python \
	--disable-gcc-flags \
	--enable-silent-rules
make -j ${CPU_COUNT}
make -j ${CPU_COUNT} check
make install
