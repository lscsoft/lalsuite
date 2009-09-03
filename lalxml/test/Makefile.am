AMCFLAGS = $(LALSUPPORT_CFLAGS)

if LALSUPPORT
check_PROGRAMS = LALXMLTest
else
check_PROGRAMS =
endif

pkgdata_DATA = VOTable-1.1.xsd
EXTRA_DIST = VOTable-1.1.xsd

LALXMLTest_SOURCES = LALXMLTest.c
LALXMLTest_LDADD = $(top_builddir)/src/liblalxml.la \
	$(LALSUPPORT_LIBS)

INCLUDES = -I../include

if LALSUPPORT
TESTS = $(check_PROGRAMS)
else
TESTS =
endif

TESTS_ENVIRONMENT = 1>/dev/null
