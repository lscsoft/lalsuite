## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libframedata.la

OLDFRAMEDATASOURCES = FrameData.c SpecBuffer.c DataBuffer.c
FRAMEDATASOURCES = FrameStream.c

if FRAME
libframedata_la_SOURCES = $(FRAMEDATASOURCES)
EXTRA_TEXSOURCES =
else
libframedata_la_SOURCES = 
EXTRA_TEXSOURCES = $(FRAMEDATASOURCES)
endif

EXTRA_DIST = $(FRAMEDATASOURCES) $(OLDFRAMEDATASOURCES)

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
