## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = AVFactories.h SeqFactories.h

all-local: $(top_builddir)/include/lal/AVFactories.h $(top_builddir)/include/lal/SeqFactories.h

$(top_builddir)/include/lal/AVFactories.h: AVFactories.h
	cp -f $(srcdir)/AVFactories.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/SeqFactories.h: SeqFactories.h
	cp -f $(srcdir)/SeqFactories.h $(top_builddir)/include/lal
