## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = Overlap.h Dirichlet.h CrossCorr.h SiteInfo.h

all-local: $(top_builddir)/include/lal/Overlap.h   \
           $(top_builddir)/include/lal/Dirichlet.h \
           $(top_builddir)/include/lal/CrossCorr.h \
           $(top_builddir)/include/lal/SiteInfo.h

$(top_builddir)/include/lal/Overlap.h: Overlap.h
	cp -f $(srcdir)/Overlap.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/Dirichlet.h: Dirichlet.h
	cp -f $(srcdir)/Dirichlet.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/CrossCorr.h: CrossCorr.h
	cp -f $(srcdir)/CrossCorr.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/SiteInfo.h: SiteInfo.h
	cp -f $(srcdir)/SiteInfo.h $(top_builddir)/include/lal
