## Process this file with automake to produce Makefile.in
check_PROGRAMS = TestConstructPLUT TestDriveHough TestHoughMap TestPeak2PHMD
TestConstructPLUT_SOURCES =  TestConstructPLUT.c
TestConstructPLUT_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TestDriveHough_SOURCES = TestDriveHough.c
TestDriveHough_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TestHoughMap_SOURCES = TestHoughMap.c
TestHoughMap_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TestPeak2PHMD_SOURCES = TestPeak2PHMD.c
TestPeak2PHMD_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS =  TestConstructPLUT TestDriveHough TestHoughMap TestPeak2PHMD
MOSTLYCLEANFILES = .dvi-dep OutHough.asc
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
