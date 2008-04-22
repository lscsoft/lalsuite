## Process this file with automake to produce Makefile.in
# $Id$

# for using lalapps-defines + functions
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps  -I$(top_srcdir)/src/pulsar/FDS_isolated
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_pulsar_radiometer

lalapps_pulsar_radiometer_SOURCES = radiometer_aux.c radiometer_sfts.c radiometer.h 
