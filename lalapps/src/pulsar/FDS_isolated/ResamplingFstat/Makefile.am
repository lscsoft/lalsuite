## Process this file with automake to produce Makefile.in
# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps

## make sure lalapps.la is built
lalappsla = $(top_builddir)/src/lalapps/liblalapps.la
$(lalappsla):
	cd $(top_builddir)/src/lalapps && $(MAKE)
LDADD = $(lalappsla)


bin_PROGRAMS = lalapps_ComputeFStatistic_resamp

EXTRA_PROGRAMS = resampling

resampling_SOURCES = proc.c fstat.c

lalapps_ComputeFStatistic_resamp_SOURCES = ComputeFStatistic_resamp.c ../HeapToplist.c ../HeapToplist.h
