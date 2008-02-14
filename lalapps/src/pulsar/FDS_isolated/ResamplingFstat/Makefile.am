## Process this file with automake to produce Makefile.in
# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_ComputeFStatistic_resamp

EXTRA_PROGRAMS = resampling

resampling_SOURCES = proc.c fstat.c

lalapps_ComputeFStatistic_resamp_SOURCES = ComputeFStatistic_resamp.c ../HeapToplist.c ../HeapToplist.h
