## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = Random.h FindRoot.h Integrate.h Interpolate.h Sort.h

all-local: $(top_builddir)/include/lal/Random.h   \
           $(top_builddir)/include/lal/FindRoot.h \
           $(top_builddir)/include/lal/Integrate.h \
           $(top_builddir)/include/lal/Interpolate.h \
           $(top_builddir)/include/lal/Sort.h

$(top_builddir)/include/lal/Random.h: Random.h
	cp -f $(srcdir)/Random.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/FindRoot.h: FindRoot.h
	cp -f $(srcdir)/FindRoot.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/Integrate.h: Integrate.h
	cp -f $(srcdir)/Integrate.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/Interpolate.h: Interpolate.h
	cp -f $(srcdir)/Interpolate.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/Sort.h: Sort.h
	cp -f $(srcdir)/Sort.h $(top_builddir)/include/lal
