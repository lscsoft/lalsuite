## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libframedata.la

OLD_FRAMEDATA_SOURCES = FrameData.c SpecBuffer.c DataBuffer.c
FRAMEDATA_SOURCES = FrameStream.c

if FRAME
libframedata_la_SOURCES = $(FRAMEDATA_SOURCES)
EXTRA_TEXSOURCES =
else
libframedata_la_SOURCES = 
EXTRA_TEXSOURCES = $(FRAMEDATA_SOURCES)
endif

EXTRA_DIST = $(FRAMEDATA_SOURCES) $(OLD_FRAMEDATA_SOURCES)

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
