## Process this file with automake to produce Makefile.in
check_PROGRAMS = BandPassTest IIRFilterTest LPCTest
BandPassTest_SOURCES = BandPassTest.c
BandPassTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
IIRFilterTest_SOURCES = IIRFilterTest.c
IIRFilterTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LPCTest_SOURCES = LPCTest.c
LPCTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = BandPassTest IIRFilterTest LPCTest
TESTS_ENVIRONMENT = 1>/dev/null
MOSTLYCLEANFILES = .dvi-dep out.dat
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
