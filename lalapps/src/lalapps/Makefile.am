## Process this file with automake to produce Makefile.in
# $Id$
if FRAME
FRAMESRCS = series.c series.h lalappsfrutils.c lalappsfrutils.h
else
FRAMESRCS = 
endif
pythonlibdir = ${prefix}/lib/python
pythonlib_SCRIPTS = pipeline.py
noinst_LTLIBRARIES = liblalapps.la
liblalapps_la_SOURCES = lalapps.c lalapps.h getopt.c getopt1.c getopt.h \
                        processtable.c processtable.h $(FRAMESRCS)
INCLUDES = -I$(top_srcdir)/src
EXTRA_DIST = series.c series.h $(SCRIPTS)
