## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = Random.h FindRoot.h Integrate.h Interpolate.h Sort.h

all-local: $(top_srcdir)/include/Random.h   \
           $(top_srcdir)/include/FindRoot.h \
           $(top_srcdir)/include/Integrate.h \
           $(top_srcdir)/include/Interpolate.h \
           $(top_srcdir)/include/Sort.h

$(top_srcdir)/include/Random.h: Random.h
	cp -f Random.h $(top_srcdir)/include

$(top_srcdir)/include/FindRoot.h: FindRoot.h
	cp -f FindRoot.h $(top_srcdir)/include

$(top_srcdir)/include/Integrate.h: Integrate.h
	cp -f Integrate.h $(top_srcdir)/include

$(top_srcdir)/include/Interpolate.h: Interpolate.h
	cp -f Interpolate.h $(top_srcdir)/include

$(top_srcdir)/include/Sort.h: Sort.h
	cp -f Sort.h $(top_srcdir)/include
