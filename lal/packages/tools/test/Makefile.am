## Process this file with automake to produce Makefile.in
check_PROGRAMS = UnitsTest DetectorSiteTest LALTestDetResponse0
UnitsTest_SOURCES = UnitsTest.c
UnitsTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
DetectorSiteTest_SOURCES = DetectorSiteTest.c
DetectorSiteTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALTestDetResponse0_SOURCES = LALTestDetResponse0.c
LALTestDetResponse0_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = UnitsTest DetectorSiteTest LALTestDetResponse0
MOSTLYCLEANFILES = .dvi-dep PrintVector.*
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
