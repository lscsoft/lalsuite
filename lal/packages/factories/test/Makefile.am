## Process this file with automake to produce Makefile.in
check_PROGRAMS = \
  VectorFactoriesTest \
  VectorSequenceFactoriesTest \
  ArrayFactoriesTest

VectorFactoriesTest_SOURCES = VectorFactoriesTest.c
VectorSequenceFactoriesTest_SOURCES = VectorSequenceFactoriesTest.c

VectorFactoriesTest_LDADD = $(top_srcdir)/lib/liblal.la
VectorSequenceFactoriesTest_LDADD = $(top_srcdir)/lib/liblal.la

ArrayFactoriesTest_SOURCES = ArrayFactoriesTest.c
ArrayFactoriesTest_LDADD = $(top_srcdir)/lib/liblal.la

INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include
TESTS = VectorFactoriesTest \
        VectorSequenceFactoriesTest \
        ArrayFactoriesTest

EXTRA_DIST = \
        VectorFactoriesTest.m4 VectorFactoriesTestFunction.m4 \
        VectorSequenceFactoriesTest.m4 VectorSequenceFactoriesTestFunction.m4 \
        ArrayFactoriesTest.m4 ArrayFactoriesTestFunction.m4
MAINTAINERCLEANFILES = $(SOURCES)

VectorFactoriesTest.c: VectorFactoriesTest.m4 VectorFactoriesTestFunction.m4
	m4 VectorFactoriesTest.m4 > VectorFactoriesTest.c

VectorSequenceFactoriesTest.c: VectorSequenceFactoriesTest.m4 VectorSequenceFactoriesTestFunction.m4
	m4 VectorSequenceFactoriesTest.m4 > VectorSequenceFactoriesTest.c

ArrayFactoriesTest.c: ArrayFactoriesTest.m4 ArrayFactoriesTestFunction.m4
	m4 ArrayFactoriesTest.m4 > ArrayFactoriesTest.c
