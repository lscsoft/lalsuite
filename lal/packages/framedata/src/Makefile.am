## Process this file with automake to produce Makefile.in
OLDFRAMEDATASOURCES = FrameData.c SpecBuffer.c DataBuffer.c
FRAMEDATASOURCES = FrameStream.c FrameCache.c FrameStreamDef.h \
                   FrameCalibration.c

if FRAME
lib_LTLIBRARIES = liblalframe.la
liblalframe_la_SOURCES = $(FRAMEDATASOURCES) FrameSeries.c
liblalframe_la_LIBADD = $(FRAME_LIBS)
EXTRA_TEXSOURCES =
else
EXTRA_TEXSOURCES = $(FRAMEDATASOURCES) FrameSeries.c
endif

EXTRA_DIST = $(FRAMEDATASOURCES) $(OLDFRAMEDATASOURCES) \
        FrameSeries.c FrameSeries.m4 FrameSeriesRead.m4 FrameSeriesWrite.m4

INCLUDES = $(FRAME_CFLAGS) -I$(top_builddir)/include -I$(srcdir)/../include
MAINTAINERCLEANFILES = FrameSeries.c
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep

FrameSeries.c: FrameSeries.m4 FrameSeriesRead.m4 FrameSeriesWrite.m4
	cd $(srcdir) && m4 FrameSeries.m4 > FrameSeries.c
@DVIMKFILE@
