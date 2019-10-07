#!/bin/bash

# when running on gitlab-ci, we are not using a production
# build, so we don't want to use NDEBUG
export CPPFLAGS="${CPPFLAGS} -UNDEBUG"

./configure \
	--prefix=${PREFIX} \
	--enable-help2man \
	--enable-cfitsio \
	--enable-openmp \
	--enable-mpi \
	MPICC=${PREFIX}/bin/mpicc \
	MPICXX=${PREFIX}/bin/mpicxx \
	MPIFC=${PREFIX}/bin/mpifc

make -j ${CPU_COUNT}
make -j ${CPU_COUNT} check
make -j ${CPU_COUNT} install
