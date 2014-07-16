include $(top_srcdir)/gnuscripts/lalsuite_header_links.am

TEST_CPPFLAGS = -DTEST_DATA_DIR='"$(abs_srcdir)/"'

check_PROGRAMS = \
	NoisePSDTest

TESTS = \
	$(check_PROGRAMS)

EXTRA_DIST = NoisePSDTest.eps
EXTRA_TEXSOURCES = NoisePSDTest.eps

CLEANFILES = *.out
