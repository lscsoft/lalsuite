## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = LALfct.h LALFCTInterface.h

all-local: $(top_builddir)/include/lal/LALfct.h

$(top_builddir)/include/lal/LALfct.h: LALfct.h
	cp -f $(srcdir)/LALfct.h $(top_builddir)/include/lal

all-local: $(top_builddir)/include/lal/LALFCTInterface.h

$(top_builddir)/include/lal/LALFCTInterface.h: LALFCTInterface.h
	cp -f $(srcdir)/LALFCTInterface.h $(top_builddir)/include/lal
