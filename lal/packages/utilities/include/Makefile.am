## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = Random.h FindRoot.h

all-local: $(top_srcdir)/include/Random.h   \
           $(top_srcdir)/include/FindRoot.h

$(top_srcdir)/include/Random.h: Random.h
	cp -f Random.h $(top_srcdir)/include

$(top_srcdir)/include/FindRoot.h: FindRoot.h
	cp -f FindRoot.h $(top_srcdir)/include
