## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = BandPassTimeSeries.h IIRFilter.h ZPGFilter.h

all-local: $(top_builddir)/include/lal/BandPassTimeSeries.h   \
           $(top_builddir)/include/lal/IIRFilter.h \
           $(top_builddir)/include/lal/ZPGFilter.h

$(top_builddir)/include/lal/BandPassTimeSeries.h: BandPassTimeSeries.h
	cp -f $(srcdir)/BandPassTimeSeries.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/IIRFilter.h: IIRFilter.h
	cp -f $(srcdir)/IIRFilter.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/ZPGFilter.h: ZPGFilter.h
	cp -f $(srcdir)/ZPGFilter.h $(top_builddir)/include/lal
