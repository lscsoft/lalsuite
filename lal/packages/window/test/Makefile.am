include $(top_srcdir)/gnuscripts/lallibs_header_links.am

AM_CPPFLAGS = -I$(top_builddir)/include -I$(srcdir)/../include @AM_CPPFLAGS@
LDADD = $(top_builddir)/lib/lal/liblal.la -lm

check_PROGRAMS = WindowTest

TESTS = $(check_PROGRAMS)
