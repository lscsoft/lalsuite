## Process this file with automake to produce Makefile.in
check_PROGRAMS = SlopeDetectorFilterTest
SlopeDetectorFilterTest_SOURCES = SlopeDetectorFilterTest.c
SlopeDetectorFilterTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = SlopeDetectorFilterTest
TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
