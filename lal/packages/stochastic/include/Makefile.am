## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = Overlap.h Dirichlet.h

all-local: $(top_srcdir)/include/Overlap.h   \
           $(top_srcdir)/include/Dirichlet.h

$(top_srcdir)/include/Overlap.h: Overlap.h
	cp -f Overlap.h $(top_srcdir)/include

$(top_srcdir)/include/Dirichlet.h: Dirichlet.h
	cp -f Dirichlet.h $(top_srcdir)/include
