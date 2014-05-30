include $(top_srcdir)/gnuscripts/lallibs_header_links.am

TEST_CPPFLAGS = -DTEST_DATA_DIR='"$(abs_srcdir)/"'
LDADD = $(top_builddir)/lib/lal/liblal.la -lm

check_PROGRAMS = WindowTest

TESTS = $(check_PROGRAMS)
