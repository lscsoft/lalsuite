## Process this file with automake to produce Makefile.in
# $Id$
bin_PROGRAMS = lalapps_detresponse
lalapps_detresponse_SOURCES = detresponse.c util.c set_source_params.c generate_timeseries_response.c cmdline.c getopt.c getopt1.c mystrncasecmp.c
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I.
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
