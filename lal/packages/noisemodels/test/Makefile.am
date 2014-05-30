include $(top_srcdir)/gnuscripts/lallibs_header_links.am

TEST_CPPFLAGS = -DTEST_DATA_DIR='"$(abs_srcdir)/"'
LDADD = \
	$(top_builddir)/lib/lal/liblal.la \
	$(top_builddir)/lib/lalsupport/src/liblalsupport.la

check_PROGRAMS = \
	NoisePSDTest

TESTS = \
	$(check_PROGRAMS)

EXTRA_DIST = NoisePSDTest.eps
EXTRA_TEXSOURCES = NoisePSDTest.eps

CLEANFILES = *.out
