INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

noinst_HEADERS = \
	ArrayFactoriesTest_source.c \
	ArraySequenceFactoriesTest_source.c \
	VectorFactoriesTest_source.c \
	VectorSequenceFactoriesTest_source.c

check_PROGRAMS = \
	ArrayFactoriesTest \
	ArraySequenceFactoriesTest \
	VectorFactoriesTest \
	StringVectorTest \
	VectorSequenceFactoriesTest

TESTS = \
	$(check_PROGRAMS)
