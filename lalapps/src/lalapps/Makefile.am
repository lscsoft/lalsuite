## Process this file with automake to produce Makefile.in
# $Id$
if FRAME
FRAMESRCS = series.c series.h lalappsfrutils.c lalappsfrutils.h
else
FRAMESRCS = 
endif
noinst_LTLIBRARIES = liblalapps.la
liblalapps_la_SOURCES = lalapps.c lalapps.h getopt.c getopt1.c getopt.h \
                        processtable.c processtable.h $(FRAMESRCS)
bin_PROGRAMS = lalapps_version
lalapps_version_SOURCES = version.c lalapps.h
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
INCLUDES = -I$(top_srcdir)/src
EXTRA_DIST = series.c series.h
TESTS = lalapps_version
