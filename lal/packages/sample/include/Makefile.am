## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = LALSample.h LALTemplate.h

all-local: $(top_builddir)/include/LALSample.h   \
           $(top_builddir)/include/LALTemplate.h

$(top_builddir)/include/LALSample.h: LALSample.h
	cp -f $(srcdir)/LALSample.h $(top_builddir)/include

$(top_builddir)/include/LALTemplate.h: LALTemplate.h
	cp -f $(srcdir)/LALTemplate.h $(top_builddir)/include
