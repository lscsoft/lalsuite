## Process this file with automake to produce Makefile.in
check_PROGRAMS = SlopeDetectorFilterTest
SlopeDetectorFilterTest_SOURCES = SlopeDetectorFilterTest.c
SlopeDetectorFilterTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = SlopeDetectorFilterTest
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
