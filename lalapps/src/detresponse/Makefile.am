## Process this file with automake to produce Makefile.in
# $Id$
bin_PROGRAMS = lalapps_detresponse
lalapps_detresponse_SOURCES = detresponse.c detresponse.h util.c util.h set_source_params.c generate_timeseries_response.c cmdline.c cmdline.h
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I.
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
