## Process this file with automake to produce Makefile.in
check_PROGRAMS = TestConstructPLUT TestDriveHough TestHoughMap TestPeak2PHMD \
	TestNDConstructPLUT TestDriveNDHough TestNDHoughMap TestNDPeak2PHMD \
	TestStatistics TestVelocity

TestConstructPLUT_SOURCES =  TestConstructPLUT.c
TestConstructPLUT_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TestDriveHough_SOURCES = TestDriveHough.c
TestDriveHough_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TestHoughMap_SOURCES = TestHoughMap.c
TestHoughMap_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TestPeak2PHMD_SOURCES = TestPeak2PHMD.c
TestPeak2PHMD_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

TestNDConstructPLUT_SOURCES =  TestNDConstructPLUT.c
TestNDConstructPLUT_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TestDriveNDHough_SOURCES = TestDriveNDHough.c
TestDriveNDHough_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TestNDHoughMap_SOURCES = TestNDHoughMap.c
TestNDHoughMap_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TestNDPeak2PHMD_SOURCES = TestNDPeak2PHMD.c
TestNDPeak2PHMD_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

TestStatistics_SOURCES = TestStatistics.c
TestStatistics_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TestVelocity_SOURCES = TestVelocity.c
TestVelocity_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS =  TestConstructPLUT TestDriveHough TestHoughMap TestPeak2PHMD \
	TestNDConstructPLUT TestDriveNDHough TestNDHoughMap TestNDPeak2PHMD \
	TestStatistics TestVelocity


TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null LAL_DATA_PATH=$(srcdir)/../../pulsar/test
#TESTS_ENVIRONMENT = LAL_DATA_PATH=$(srcdir)
#noinst_DATA = ../../pulsar/test/earth03.dat ../../pulsar/test/sun03.dat
#EXTRA_DIST = $(DATA)

MOSTLYCLEANFILES = .dvi-dep OutHough.asc OutHistogram.asc
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
