## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = LALVersion.h LALMalloc.h LALError.h \
                     PrintVector.h PrintFTSeries.h

all-local: $(top_builddir)/include/lal/LALVersion.h \
           $(top_builddir)/include/lal/LALMalloc.h \
           $(top_builddir)/include/lal/LALError.h \
           $(top_builddir)/include/lal/PrintVector.h \
           $(top_builddir)/include/lal/PrintFTSeries.h

$(top_builddir)/include/lal/LALVersion.h: LALVersion.h
	cp -f $(srcdir)/LALVersion.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/LALMalloc.h: LALMalloc.h
	cp -f $(srcdir)/LALMalloc.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/LALError.h: LALError.h
	cp -f $(srcdir)/LALError.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/PrintVector.h: PrintVector.h
	cp -f $(srcdir)/PrintVector.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/PrintFTSeries.h: PrintFTSeries.h
	cp -f $(srcdir)/PrintFTSeries.h $(top_builddir)/include/lal
