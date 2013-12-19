AM_CFLAGS = $(LALSUPPORT_CFLAGS) -DTEST_DATA_DIR='"$(abs_srcdir)/"' @AM_CFLAGS@
AM_CPPFLAGS = -I$(top_builddir)/include @AM_CPPFLAGS@
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalxml.la

check_PROGRAMS = LALXMLTest
TESTS = $(check_PROGRAMS)

pkgdata_DATA = VOTable-1.1.xsd
EXTRA_DIST = VOTable-1.1.xsd
