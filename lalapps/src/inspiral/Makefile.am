## Process this file with automake to produce Makefile.in
# $Id$
INSPIRALSRC = inspiral.c inspiral.h
if FRAME
bin_PROGRAMS = lalapps_inspiral
lalapps_inspiral_SOURCES = $(INSPIRALSRC)
else
bin_PROGRAMS = 
endif
EXTRA_DIST = $(INSPIRALSRC)

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

MOSTLYCLEANFILES = results.xml
