## Process this file with automake to produce Makefile.in
# $Id$
if FRAME
FRAMESRCS = series.c series.h
else
FRAMESRCS = 
endif
noinst_LTLIBRARIES = liblalapps.la
liblalapps_la_SOURCES = lalapps.c lalapps.h \
                        processtable.c processtable.h $(FRAMESRCS)
INCLUDES = -I$(top_srcdir)/src
EXTRA_DIST = series.c series.h
