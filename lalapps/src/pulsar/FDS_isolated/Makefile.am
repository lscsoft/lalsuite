## Process this file with automake to produce Makefile.in
# $Id$

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

## build everything completely static
AM_LDFLAGS = -all-static

if FRAME
FRAMEPROGRAMS = lalapps_FrComputeFStatistic
else
FRAMEPROGRAMS = 
endif

bin_PROGRAMS = lalapps_ComputeFStatistic lalapps_FstatShapeTest lalapps_SemiAnalyticF lalapps_polka lalapps_FstatShapeTestLAL ComputeFStatisticBinary  compareFstats $(FRAMEPROGRAMS) ComputeFStatistic_v2

lalapps_ComputeFStatistic_SOURCES = ComputeFStatistic.c ComputeFStatistic.h	\
	rngmed.c rngmed.h clusters.c clusters.h  			\
	DopplerScan.c DopplerScan.h NeighboringGridPoints.c

lalapps_FrComputeFStatistic_SOURCES = FrComputeFStatistic.c FrComputeFStatistic.h	\
	rngmed.c rngmed.h clusters.c clusters.h  			\
	DopplerScan.c DopplerScan.h NeighboringGridPoints.c

lalapps_FstatShapeTest_SOURCES = FstatShapeTest.c 
lalapps_FstatShapeTestLAL_SOURCES = FstatShapeTestLAL.c

lalapps_SemiAnalyticF_SOURCES =	SemiAnalyticF.c SemiAnalyticF.h

lalapps_polka_SOURCES = polka.c

ComputeFStatisticBinary_SOURCES = ComputeFStatisticBinary.c ComputeFStatisticBinary.h rngmed.c rngmed.h clusters.c clusters.h

compareFstats_SOURCES = compareFstats.c

ComputeFStatistic_v2_SOURCES = ComputeFStatistic_v2.c 	rngmed.c rngmed.h clusters.c clusters.h  			\
	DopplerScan.c DopplerScan.h NeighboringGridPoints.c