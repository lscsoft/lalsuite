## Process this file with automake to produce Makefile.in
# $Id$
bin_PROGRAMS = lalapps_inspiral

lalapps_inspiral_SOURCES = inspiral.c inspiral.h

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

MOSTLYCLEANFILES = results.xml
