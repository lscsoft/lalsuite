## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = BandPassTimeSeries.h IIRFilter.h ZPGFilter.h

all-local: $(top_builddir)/include/BandPassTimeSeries.h   \
           $(top_builddir)/include/IIRFilter.h \
           $(top_builddir)/include/ZPGFilter.h

$(top_builddir)/include/BandPassTimeSeries.h: BandPassTimeSeries.h
	cp -f $(srcdir)/BandPassTimeSeries.h $(top_builddir)/include

$(top_builddir)/include/IIRFilter.h: IIRFilter.h
	cp -f $(srcdir)/IIRFilter.h $(top_builddir)/include

$(top_builddir)/include/ZPGFilter.h: ZPGFilter.h
	cp -f $(srcdir)/ZPGFilter.h $(top_builddir)/include
