## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = FrameData.h SpecBuffer.h DataBuffer.h

all-local: \
        $(top_builddir)/include/lal/FrameData.h \
        $(top_builddir)/include/lal/SpecBuffer.h \
        $(top_builddir)/include/lal/DataBuffer.h

$(top_builddir)/include/lal/FrameData.h: FrameData.h
	cp -f $(srcdir)/FrameData.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/SpecBuffer.h: SpecBuffer.h
	cp -f $(srcdir)/SpecBuffer.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/DataBuffer.h: DataBuffer.h
	cp -f $(srcdir)/DataBuffer.h $(top_builddir)/include/lal
