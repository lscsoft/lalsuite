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
bin_PROGRAMS = lalapps_version lalapps_tconvert lalapps_cache
lalapps_version_SOURCES = version.c lalapps.h
lalapps_tconvert_SOURCES = tconvert.c getdate.c getdate.h
lalapps_cache_SOURCES = cache.c
LDADD = liblalapps.la @LAL_LIBS@
INCLUDES = -I$(top_srcdir)/src
if HAVE_PYTHON
pkgpython_DATA = __init__.py defaults.py
endif
EXTRA_DIST = series.c series.h getdate.y __init__.py
TESTS = lalapps_version
