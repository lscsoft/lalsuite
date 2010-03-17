LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

check_PROGRAMS = \
	GeocentricGeodeticTest \
	SkyCoordinatesTest

SkyCoordinatesTest_SOURCES = SkyCoordinatesTest.c

GeocentricGeodeticTest_SOURCES = GeocentricGeodeticTest.c

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

TESTS = \
	GeocentricGeodeticTest \
	SkyCoordinatesTest

TESTS_ENVIRONMENT = 1>/dev/null

EXTRA_DIST = \
	CheckStatus.c \
	CheckStatus.h \
	InjectTest.c

MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex

dvi-local: .dvi-dep

@DVIMKFILE@
