AM_CFLAGS = $(LALSUPPORT_CFLAGS) -DDATADIR='"$(abs_srcdir)/"'
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
