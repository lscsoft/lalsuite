## Process this file with automake to produce Makefile.in
check_PROGRAMS = UnitsTest DetectorSiteTest
UnitsTest_SOURCES = UnitsTest.c
UnitsTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
DetectorSiteTest_SOURCES = DetectorSiteTest.c
DetectorSiteTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = UnitsTest DetectorSiteTest
MOSTLYCLEANFILES = .dvi-dep PrintVector.*
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
