AM_CFLAGS = $(LALSUPPORT_CFLAGS)
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalxml.la

if LALSUPPORT
check_PROGRAMS = LALXMLTest
TESTS = $(check_PROGRAMS)
else
check_PROGRAMS =
TESTS =
endif

pkgdata_DATA = VOTable-1.1.xsd
EXTRA_DIST = VOTable-1.1.xsd

TESTS_ENVIRONMENT = 1>/dev/null LAL_DATA_PATH=$(srcdir)
