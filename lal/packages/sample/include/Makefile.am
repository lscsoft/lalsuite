## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = LALSample.h LALTemplate.h

all-local: $(top_builddir)/include/lal/LALSample.h   \
           $(top_builddir)/include/lal/LALTemplate.h

$(top_builddir)/include/lal/LALSample.h: LALSample.h
	cp -f $(srcdir)/LALSample.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/LALTemplate.h: LALTemplate.h
	cp -f $(srcdir)/LALTemplate.h $(top_builddir)/include/lal
