## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = AVFactories.h SeqFactories.h

all-local: $(top_srcdir)/include/AVFactories.h $(top_srcdir)/include/SeqFactories.h

$(top_srcdir)/include/AVFactories.h: AVFactories.h
	cp -f AVFactories.h $(top_srcdir)/include

$(top_srcdir)/include/SeqFactories.h: SeqFactories.h
	cp -f SeqFactories.h $(top_srcdir)/include
