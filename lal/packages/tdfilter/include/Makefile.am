## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = BandPassTimeSeries.h IIRFilter.h ZPGFilter.h

all-local: $(top_srcdir)/include/BandPassTimeSeries.h   \
           $(top_srcdir)/include/IIRFilter.h \
           $(top_srcdir)/include/ZPGFilter.h

$(top_srcdir)/include/BandPassTimeSeries.h: BandPassTimeSeries.h
	cp -f BandPassTimeSeries.h $(top_srcdir)/include

$(top_srcdir)/include/IIRFilter.h: IIRFilter.h
	cp -f IIRFilter.h $(top_srcdir)/include

$(top_srcdir)/include/ZPGFilter.h: ZPGFilter.h
	cp -f ZPGFilter.h $(top_srcdir)/include
