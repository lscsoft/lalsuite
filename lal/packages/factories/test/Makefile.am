LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

noinst_HEADERS = \
	ArrayFactoriesTest_source.c \
	ArraySequenceFactoriesTest_source.c \
	VectorFactoriesTest_source.c \
	VectorSequenceFactoriesTest_source.c

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
