## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = AVFactories.h SeqFactories.h

all-local: $(top_builddir)/include/AVFactories.h $(top_builddir)/include/SeqFactories.h

$(top_builddir)/include/AVFactories.h: AVFactories.h
	cp -f $(srcdir)/AVFactories.h $(top_builddir)/include

$(top_builddir)/include/SeqFactories.h: SeqFactories.h
	cp -f $(srcdir)/SeqFactories.h $(top_builddir)/include
