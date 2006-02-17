## Process this file with automake to produce Makefile.in
# $Id$
#
# Charlie
#
TRACKSEARCHSRCS = tracksearch.c
TSDATAGENSRCS = TSDatagen.c 

bin_PROGRAMS = lalapps_tracksearch lalapps_TSDatagen
lalapps_tracksearch_SOURCES = $(TRACKSEARCHSRCS)
lalapps_TSDatagen_SOURCES = $(TSDATAGENSRCS)


EXTRA_DIST = $(TRACKSEARCHSRCS)
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la





