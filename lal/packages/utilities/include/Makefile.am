## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = Random.h FindRoot.h Integrate.h Interpolate.h Sort.h

all-local: $(top_builddir)/include/Random.h   \
           $(top_builddir)/include/FindRoot.h \
           $(top_builddir)/include/Integrate.h \
           $(top_builddir)/include/Interpolate.h \
           $(top_builddir)/include/Sort.h

$(top_builddir)/include/Random.h: Random.h
	cp -f $(srcdir)/Random.h $(top_builddir)/include

$(top_builddir)/include/FindRoot.h: FindRoot.h
	cp -f $(srcdir)/FindRoot.h $(top_builddir)/include

$(top_builddir)/include/Integrate.h: Integrate.h
	cp -f $(srcdir)/Integrate.h $(top_builddir)/include

$(top_builddir)/include/Interpolate.h: Interpolate.h
	cp -f $(srcdir)/Interpolate.h $(top_builddir)/include

$(top_builddir)/include/Sort.h: Sort.h
	cp -f $(srcdir)/Sort.h $(top_builddir)/include
