## Process this file with automake to produce Makefile.in
# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
AM_LDFLAGS = $(top_builddir)/src/lalapps/liblalapps.la

if FRAME
FRAMEPROGRAMS = lalapps_FrComputeFStatistic
else
FRAMEPROGRAMS = 
endif

bin_PROGRAMS = lalapps_ComputeFStatistic lalapps_FstatShapeTest lalapps_SemiAnalyticF \
	lalapps_polka lalapps_FstatShapeTestLAL compareFstats ComputeFStatistic_v2 \
	lalapps_makeInvetofile lalapps_extractSFTband lalapps_FindSh

CFSsources = 	ComputeFStatistic.h	\
		clusters.c clusters.h  	\
		DopplerScan.c DopplerScan.h NeighboringGridPoints.c

lalapps_ComputeFStatistic_SOURCES = ComputeFStatistic.c $(CFSsources)
lalapps_ComputeFStatistic_LDFLAGS = -all-static

lalapps_FrComputeFStatistic_SOURCES = FrComputeFStatistic.c FrComputeFStatistic.h	\
	rngmed.c rngmed.h clusters.c clusters.h  			\
	DopplerScan.c DopplerScan.h NeighboringGridPoints.c

lalapps_FstatShapeTest_SOURCES = FstatShapeTest.c 
lalapps_FstatShapeTestLAL_SOURCES = FstatShapeTestLAL.c

lalapps_SemiAnalyticF_SOURCES =	SemiAnalyticF.c SemiAnalyticF.h

lalapps_polka_SOURCES = polka.c

lalapps_makeInvetofile_SOURCES=makeInvetofile.c

lalapps_extractSFTband_SOURCES=extractSFTband.c

lalapps_FindSh_SOURCES=FindSh.c FindSh.h clusters.c clusters.h

compareFstats_SOURCES = compareFstats.c

ComputeFStatistic_v2_SOURCES = ComputeFStatistic_v2.c clusters.c clusters.h  \
	DopplerScan.c DopplerScan.h NeighboringGridPoints.c

