## Process this file with automake to produce Makefile.in
check_PROGRAMS = BandPassTest IIRFilterTest
BandPassTest_SOURCES = BandPassTest.c
BandPassTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
IIRFilterTest_SOURCES = IIRFilterTest.c
IIRFilterTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = BandPassTest IIRFilterTest
MOSTLYCLEANFILES = .dvi-dep out.dat
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
