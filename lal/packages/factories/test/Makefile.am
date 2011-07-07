LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

noinst_HEADERS = \
	ArrayFactoriesTest_source.c \
	ArraySequenceFactoriesTest_source.c

check_PROGRAMS = \
	ArrayFactoriesTest \
	ArraySequenceFactoriesTest \
	VectorFactoriesTest \
	VectorSequenceFactoriesTest

TESTS = \
	ArrayFactoriesTest \
	ArraySequenceFactoriesTest \
	VectorFactoriesTest \
	VectorSequenceFactoriesTest

EXTRA_DIST = \
	VectorFactoriesTest.m4 \
	VectorFactoriesTestFunction.m4 \
	VectorSequenceFactoriesTest.m4 \
	VectorSequenceFactoriesTestFunction.m4

MAINTAINERCLEANFILES = \
	VectorFactoriesTest.c \
	VectorSequenceFactoriesTest.c

VectorFactoriesTest.c: VectorFactoriesTest.m4 VectorFactoriesTestFunction.m4
	cd $(srcdir) && @M4@ VectorFactoriesTest.m4 > VectorFactoriesTest.c

VectorSequenceFactoriesTest.c: VectorSequenceFactoriesTest.m4 VectorSequenceFactoriesTestFunction.m4
	cd $(srcdir) && @M4@ VectorSequenceFactoriesTest.m4 > VectorSequenceFactoriesTest.c
