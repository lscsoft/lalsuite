## Process this file with automake to produce Makefile.in
check_PROGRAMS = RingTest RingSearchTest LALRingDownTest LALWienerTest
RingTest_SOURCES = RingTest.c
RingTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
RingSearchTest_SOURCES = RingSearchTest.c
RingSearchTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALRingDownTest_SOURCES = LALRingDownTest.c CheckStatus.c CheckStatus.h
LALRingDownTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
LALWienerTest_SOURCES = LALWienerTest.c CheckStatus.c CheckStatus.h
LALWienerTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = RingTest RingSearchTest LALRingDownTest LALWienerTest
MOSTLYCLEANFILES = .dvi-dep *.out *.dat snr-*
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
