## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = LALMalloc.h LALError.h PrintVector.h PrintFTSeries.h

all-local: $(top_builddir)/include/LALMalloc.h \
           $(top_builddir)/include/LALError.h \
           $(top_builddir)/include/PrintVector.h \
           $(top_builddir)/include/PrintFTSeries.h

$(top_builddir)/include/LALMalloc.h: LALMalloc.h
	cp -f $(srcdir)/LALMalloc.h $(top_builddir)/include

$(top_builddir)/include/LALError.h: LALError.h
	cp -f $(srcdir)/LALError.h $(top_builddir)/include

$(top_builddir)/include/PrintVector.h: PrintVector.h
	cp -f $(srcdir)/PrintVector.h $(top_builddir)/include

$(top_builddir)/include/PrintFTSeries.h: PrintFTSeries.h
	cp -f $(srcdir)/PrintVector.h $(top_builddir)/include
