## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = FrameData.h SpecBuffer.h DataBuffer.h

all-local: \
        $(top_srcdir)/include/FrameData.h \
        $(top_srcdir)/include/SpecBuffer.h \
        $(top_srcdir)/include/DataBuffer.h

$(top_srcdir)/include/FrameData.h: FrameData.h
	cp -f FrameData.h $(top_srcdir)/include

$(top_srcdir)/include/SpecBuffer.h: SpecBuffer.h
	cp -f SpecBuffer.h $(top_srcdir)/include

$(top_srcdir)/include/DataBuffer.h: DataBuffer.h
	cp -f DataBuffer.h $(top_srcdir)/include
