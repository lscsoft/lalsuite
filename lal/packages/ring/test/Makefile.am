## Process this file with automake to produce Makefile.in
check_PROGRAMS = RingTest RingSearchTest
RingTest_SOURCES = RingTest.c
RingTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
RingSearchTest_SOURCES = RingSearchTest.c
RingSearchTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = RingTest RingSearchTest
MOSTLYCLEANFILES = .dvi-dep *.out *.dat snr-*
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
