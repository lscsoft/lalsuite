# Makefile for TDS_isolated 
# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_pulsar_parameter_estimation lalapps_heterodyne_pulsar

lalapps_pulsar_parameter_estimation_SOURCES = pulsar_parameter_estimation.c

lalapps_heterodyne_pulsar_SOURCES = heterodyne_pulsar.c