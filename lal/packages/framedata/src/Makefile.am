## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libframedata.la

if FRAME
libframedata_la_SOURCES = FrameData.c SpecBuffer.c DataBuffer.c
else
libframedata_la_SOURCES = 
endif

EXTRA_DIST       = FrameData.c SpecBuffer.c DataBuffer.c

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
