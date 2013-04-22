AM_CFLAGS = $(LALSUPPORT_CFLAGS) -DDATADIR='"$(abs_srcdir)/"'
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalxml.la

check_PROGRAMS = LALXMLTest
TESTS = $(check_PROGRAMS)

pkgdata_DATA = VOTable-1.1.xsd
EXTRA_DIST = VOTable-1.1.xsd
