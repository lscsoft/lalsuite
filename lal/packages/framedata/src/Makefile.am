## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = liblal.la

if FRAME
liblal_la_SOURCES = FrameData.c SpecBuffer.c DataBuffer.c
else
liblal_la_SOURCES = 
endif

EXTRA_DIST       = FrameData.c SpecBuffer.c DataBuffer.c

INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include
