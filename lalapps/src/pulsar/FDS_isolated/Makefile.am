## Process this file with automake to produce Makefile.in
# $Id$

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_ComputeFStatistic lalapps_FstatShapeTest 

lalapps_ComputeFStatistic_SOURCES = ComputeFStatistic.c ComputeFStatistic.h	\
	rngmed.c rngmed.h clusters.c clusters.h  			\
	DopplerScan.c DopplerScan.h 




lalapps_FstatShapeTest_SOURCES = FstatShapeTest.c 

