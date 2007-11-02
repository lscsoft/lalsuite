## Process this file with automake to produce Makefile.in
# $Id$

# for using lalapps-defines + functions
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps  -I$(top_srcdir)/src/pulsar/FDS_isolated
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

## not usable before automake-1.6:
## AM_LDFLAGS = -all-static

bin_PROGRAMS = pulsar_radiometer

pulsar_radiometer_SOURCES = radiometer_sfts.c radiometer.h 