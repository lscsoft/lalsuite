## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = FrameData.h SpecBuffer.h DataBuffer.h

all-local: \
        $(top_builddir)/include/FrameData.h \
        $(top_builddir)/include/SpecBuffer.h \
        $(top_builddir)/include/DataBuffer.h

$(top_builddir)/include/FrameData.h: FrameData.h
	cp -f $(srcdir)/FrameData.h $(top_builddir)/include

$(top_builddir)/include/SpecBuffer.h: SpecBuffer.h
	cp -f $(srcdir)/SpecBuffer.h $(top_builddir)/include

$(top_builddir)/include/DataBuffer.h: DataBuffer.h
	cp -f $(srcdir)/DataBuffer.h $(top_builddir)/include
