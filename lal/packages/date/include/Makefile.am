## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = Date.h TimeDelay.h

all-local: $(top_builddir)/include/lal/Date.h \
           $(top_builddir)/include/lal/TimeDelay.h

$(top_builddir)/include/lal/Date.h: Date.h
	cp -f $(srcdir)/Date.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/TimeDelay.h: TimeDelay.h
	cp -f $(srcdir)/TimeDelay.h $(top_builddir)/include/lal
