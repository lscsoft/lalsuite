## Process this file with automake to produce Makefile.in
check_PROGRAMS = RingTest RingSearchTest LALRingDownTest LALWienerTest AvgSpecTest
RingTest_SOURCES = RingTest.c
RingTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
RingSearchTest_SOURCES = RingSearchTest.c
RingSearchTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALRingDownTest_SOURCES = LALRingDownTest.c CheckStatus.c CheckStatus.h
LALRingDownTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALWienerTest_SOURCES = LALWienerTest.c CheckStatus.c CheckStatus.h
LALWienerTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
AvgSpecTest_SOURCES = AvgSpecTest.c
AvgSpecTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = RingTest RingSearchTest LALRingDownTest LALWienerTest AvgSpecTest
TESTS_ENVIRONMENT = 1>/dev/null
MOSTLYCLEANFILES = .dvi-dep *.out *.dat snr-*
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
