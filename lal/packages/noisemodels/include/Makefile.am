## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = LALNoiseModels.h

all-local: $(top_builddir)/include/lal/LALNoiseModels.h

$(top_builddir)/include/lal/LALNoiseModels.h: LALNoiseModels.h
	cp -f $(srcdir)/LALNoiseModels.h $(top_builddir)/include/lal
