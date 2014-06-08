include $(top_srcdir)/gnuscripts/lallibs_header_links.am

AM_CPPFLAGS = -I$(top_builddir)/include -I$(srcdir)/../include @AM_CPPFLAGS@
LDADD = \
	$(top_builddir)/lib/lal/liblal.la \
	$(top_builddir)/lib/lalsupport/src/liblalsupport.la

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
