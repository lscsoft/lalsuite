## Process this file with automake to produce Makefile.in
OLDFRAMEDATASOURCES = FrameData.c SpecBuffer.c DataBuffer.c
FRAMEDATASOURCES = FrameStream.c

if FRAME
lib_LTLIBRARIES = liblalframe.la
liblalframe_la_SOURCES = $(FRAMEDATASOURCES)
liblalframe_la_LIBADD = -lFrame
EXTRA_TEXSOURCES =
else
EXTRA_TEXSOURCES = $(FRAMEDATASOURCES)
endif

EXTRA_DIST = $(FRAMEDATASOURCES) $(OLDFRAMEDATASOURCES)

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
