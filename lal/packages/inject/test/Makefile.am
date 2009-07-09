## Process this file with automake to produce Makefile.in

if METAIO
NEED_METAIO = InjectionInterfaceTest
else
NEED_METAIO =
endif

check_PROGRAMS = 	BasicInjectTest InjectTest GeneratePPNInspiralTest SimulateTaylorCWTest \
			SkyCoordinatesTest GeocentricGeodeticTest GeneratePulsarSignalTest \
			GeneratePPNAmpCorInspiralTest  $(NEED_METAIO)

InjectTest_SOURCES      = InjectTest.c
InjectTest_LDADD        = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
BasicInjectTest_SOURCES = BasicInjectTest.c
BasicInjectTest_LDADD   = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
GeneratePPNInspiralTest_SOURCES = GeneratePPNInspiralTest.c
GeneratePPNInspiralTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
SimulateTaylorCWTest_SOURCES = SimulateTaylorCWTest.c
SimulateTaylorCWTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
SkyCoordinatesTest_SOURCES = SkyCoordinatesTest.c
SkyCoordinatesTest_LDADD   = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
GeocentricGeodeticTest_SOURCES = GeocentricGeodeticTest.c
GeocentricGeodeticTest_LDADD   = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
GeneratePulsarSignalTest_SOURCES = GeneratePulsarSignalTest.c
GeneratePulsarSignalTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
GeneratePPNAmpCorInspiralTest_SOURCES = GeneratePPNAmpCorInspiralTest.c
GeneratePPNAmpCorInspiralTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

if METAIO
InjectionInterfaceTest_SOURCES = InjectionInterfaceTest.c
InjectionInterfaceTest_LDADD   = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la \
			$(top_builddir)/packages/support/src/liblalmetaio.la
endif

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

TESTS = BasicInjectTest GeneratePPNInspiralTest SimulateTaylorCWTest SkyCoordinatesTest GeocentricGeodeticTest GeneratePulsarSignalTest $(NEED_MEATIO)
TESTS_ENVIRONMENT = 1>/dev/null LAL_DATA_PATH=$(top_srcdir)/packages/pulsar/test

EXTRA_DIST = CheckStatus.c CheckStatus.h L1.dat H2.dat injection.xml
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
