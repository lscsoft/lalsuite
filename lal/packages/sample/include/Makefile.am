## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = LALSample.h LALTemplate.h

all-local: $(top_srcdir)/include/LALSample.h   \
           $(top_srcdir)/include/LALTemplate.h

$(top_srcdir)/include/LALSample.h: LALSample.h
	cp -f LALSample.h $(top_srcdir)/include

$(top_srcdir)/include/LALTemplate.h: LALTemplate.h
	cp -f LALTemplate.h $(top_srcdir)/include
