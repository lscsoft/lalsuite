LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

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
	ArrayFactoriesTest.m4 \
	ArrayFactoriesTestFunction.m4 \
	ArraySequenceFactoriesTest.m4 \
	ArraySequenceFactoriesTestFunction.m4 \
	VectorFactoriesTest.m4 \
	VectorFactoriesTestFunction.m4 \
	VectorSequenceFactoriesTest.m4 \
	VectorSequenceFactoriesTestFunction.m4

MAINTAINERCLEANFILES = \
	ArrayFactoriesTest.c \
	ArraySequenceFactoriesTest.c \
	VectorFactoriesTest.c \
	VectorSequenceFactoriesTest.c

ArrayFactoriesTest.c: ArrayFactoriesTest.m4 ArrayFactoriesTestFunction.m4
	cd $(srcdir) && @M4@ ArrayFactoriesTest.m4 > ArrayFactoriesTest.c

ArraySequenceFactoriesTest.c: ArraySequenceFactoriesTest.m4 ArraySequenceFactoriesTestFunction.m4
	cd $(srcdir) && @M4@ ArraySequenceFactoriesTest.m4 > ArraySequenceFactoriesTest.c

VectorFactoriesTest.c: VectorFactoriesTest.m4 VectorFactoriesTestFunction.m4
	cd $(srcdir) && @M4@ VectorFactoriesTest.m4 > VectorFactoriesTest.c

VectorSequenceFactoriesTest.c: VectorSequenceFactoriesTest.m4 VectorSequenceFactoriesTestFunction.m4
	cd $(srcdir) && @M4@ VectorSequenceFactoriesTest.m4 > VectorSequenceFactoriesTest.c
