## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libframedata.la

FRAMEDATA_SOURCES = FrameData.c SpecBuffer.c DataBuffer.c

if FRAME
libframedata_la_SOURCES = $(FRAMEDATA_SOURCES)
else
libframedata_la_SOURCES = 
endif

EXTRA_DIST = $(FRAMEDATA_SOURCES)
EXTRA_TEXSOURCES = $(FRAMEDATA_SOURCES)

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
