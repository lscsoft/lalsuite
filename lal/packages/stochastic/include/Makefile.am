## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = Overlap.h Dirichlet.h CrossCorr.h SiteInfo.h

all-local: $(top_builddir)/include/Overlap.h   \
           $(top_builddir)/include/Dirichlet.h \
           $(top_builddir)/include/CrossCorr.h \
           $(top_builddir)/include/SiteInfo.h

$(top_builddir)/include/Overlap.h: Overlap.h
	cp -f $(srcdir)/Overlap.h $(top_builddir)/include

$(top_builddir)/include/Dirichlet.h: Dirichlet.h
	cp -f $(srcdir)/Dirichlet.h $(top_builddir)/include

$(top_builddir)/include/CrossCorr.h: CrossCorr.h
	cp -f $(srcdir)/CrossCorr.h $(top_builddir)/include

$(top_builddir)/include/SiteInfo.h: SiteInfo.h
	cp -f $(srcdir)/SiteInfo.h $(top_builddir)/include
