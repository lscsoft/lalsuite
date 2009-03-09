## Process this file with automake to produce Makefile.in
# $Id$

# for using lalapps-defines + functions
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps  -I$(top_srcdir)/src/pulsar/FDS_isolated
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_pulsar_crosscorr 

lalapps_pulsar_crosscorr_SOURCES = pulsar_crosscorr.c pulsar_crosscorr.h

