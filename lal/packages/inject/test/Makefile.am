check_PROGRAMS = \
	BasicInjectTest \
	GeneratePPNAmpCorInspiralTest \
	GeneratePPNInspiralTest \
	GeocentricGeodeticTest \
	SkyCoordinatesTest

BasicInjectTest_SOURCES = BasicInjectTest.c
BasicInjectTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

GeneratePPNInspiralTest_SOURCES = GeneratePPNInspiralTest.c
GeneratePPNInspiralTest_LDADD = \
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

GeneratePPNAmpCorInspiralTest_SOURCES = GeneratePPNAmpCorInspiralTest.c
GeneratePPNAmpCorInspiralTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

TESTS = \
	BasicInjectTest \
	GeneratePPNInspiralTest \
	GeocentricGeodeticTest \
	SkyCoordinatesTest

TESTS_ENVIRONMENT = 1>/dev/null

EXTRA_DIST = \
	CheckStatus.c \
	CheckStatus.h \
	InjectionInterfaceTest.c \
	InjectTest.c \
	injection.xml

MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex

dvi-local: .dvi-dep

@DVIMKFILE@
