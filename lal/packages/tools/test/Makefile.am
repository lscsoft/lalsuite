## Process this file with automake to produce Makefile.in
check_PROGRAMS = UnitsTest DetectorSiteTest LALTestDetResponse0 ComputeTransferTest
UnitsTest_SOURCES = UnitsTest.c
UnitsTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
DetectorSiteTest_SOURCES = DetectorSiteTest.c
DetectorSiteTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALTestDetResponse0_SOURCES = LALTestDetResponse0.c
LALTestDetResponse0_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
ComputeTransferTest_SOURCES = ComputeTransferTest.c
ComputeTransferTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = UnitsTest DetectorSiteTest LALTestDetResponse0 ComputeTransferTest
MOSTLYCLEANFILES = .dvi-dep PrintVector.* transfer.out
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
