## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = Units.h DetectorSite.h

all-local: $(top_builddir)/include/lal/Units.h \
           $(top_builddir)/include/lal/DetectorSite.h

$(top_builddir)/include/lal/Units.h: Units.h
	cp -f $(srcdir)/Units.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/DetectorSite.h: DetectorSite.h
	cp -f $(srcdir)/DetectorSite.h $(top_builddir)/include/lal
