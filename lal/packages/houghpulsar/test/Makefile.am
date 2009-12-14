INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

check_PROGRAMS = \
	TestConstructPLUT \
	TestDriveHough \
	TestDriveNDHough \
	TestHoughMap \
	TestNDConstructPLUT \
	TestNDHoughMap \
	TestNDPeak2PHMD \
	TestPeak2PHMD \
	TestStatistics \
	TestVelocity

TestConstructPLUT_SOURCES = TestConstructPLUT.c
TestConstructPLUT_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

TestDriveHough_SOURCES = TestDriveHough.c
TestDriveHough_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

TestDriveNDHough_SOURCES = TestDriveNDHough.c
TestDriveNDHough_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

TestHoughMap_SOURCES = TestHoughMap.c
TestHoughMap_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

TestNDConstructPLUT_SOURCES = TestNDConstructPLUT.c
TestNDConstructPLUT_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

TestNDHoughMap_SOURCES = TestNDHoughMap.c
TestNDHoughMap_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

TestNDPeak2PHMD_SOURCES = TestNDPeak2PHMD.c
TestNDPeak2PHMD_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

TestPeak2PHMD_SOURCES = TestPeak2PHMD.c
TestPeak2PHMD_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

TestStatistics_SOURCES = TestStatistics.c
TestStatistics_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

TestVelocity_SOURCES = TestVelocity.c
TestVelocity_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

TESTS = \
	TestConstructPLUT \
	TestDriveHough \
	TestDriveNDHough \
	TestHoughMap \
	TestNDConstructPLUT \
	TestNDHoughMap \
	TestNDPeak2PHMD \
	TestPeak2PHMD \
	TestStatistics \
	TestVelocity

TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null LAL_DATA_PATH=$(srcdir)/../../pulsar/test

DISTCLEANFILES = *.tex
MOSTLYCLEANFILES = \
	.dvi-dep \
	OutHistogram.asc \
	OutHough.asc

dvi-local: .dvi-dep

@DVIMKFILE@
