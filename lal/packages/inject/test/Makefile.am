check_PROGRAMS = \
	BasicInjectTest \
	GeneratePPNAmpCorInspiralTest \
	GeneratePPNInspiralTest \
	GeneratePulsarSignalTest \
	GeocentricGeodeticTest \
	InjectTest \
	SimulateTaylorCWTest \
	SkyCoordinatesTest

InjectTest_SOURCES = InjectTest.c
InjectTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

BasicInjectTest_SOURCES = BasicInjectTest.c
BasicInjectTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

GeneratePPNInspiralTest_SOURCES = GeneratePPNInspiralTest.c
GeneratePPNInspiralTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

SimulateTaylorCWTest_SOURCES = SimulateTaylorCWTest.c
SimulateTaylorCWTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

SkyCoordinatesTest_SOURCES = SkyCoordinatesTest.c
SkyCoordinatesTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

GeocentricGeodeticTest_SOURCES = GeocentricGeodeticTest.c
GeocentricGeodeticTest_LDADD  = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

GeneratePulsarSignalTest_SOURCES = GeneratePulsarSignalTest.c
GeneratePulsarSignalTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

GeneratePPNAmpCorInspiralTest_SOURCES = GeneratePPNAmpCorInspiralTest.c
GeneratePPNAmpCorInspiralTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

TESTS = \
	BasicInjectTest \
	GeneratePPNInspiralTest \
	GeneratePulsarSignalTest \
	GeocentricGeodeticTest \
	SimulateTaylorCWTest \
	SkyCoordinatesTest \
	$(NEED_MEATIO)

TESTS_ENVIRONMENT = 1>/dev/null LAL_DATA_PATH=$(top_srcdir)/packages/pulsar/test

EXTRA_DIST = \
	CheckStatus.c \
	CheckStatus.h \
	InjectionInterfaceTest.c \
	injection.xml

MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex

dvi-local: .dvi-dep

@DVIMKFILE@
