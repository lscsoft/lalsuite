## Process this file with automake to produce Makefile.in
# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la 

## not useable before automake-1.6
#AM_LDFLAGS = -all-static


if FRAME
FRAMEPROGRAMS = lalapps_FrComputeFStatistic
else
FRAMEPROGRAMS = 
endif

bin_PROGRAMS = lalapps_ComputeFStatistic lalapps_FstatShapeTest lalapps_SemiAnalyticF \
	lalapps_polka lalapps_uberpolka lalapps_FstatShapeTestLAL compareFstats\
	ComputeFStatistic_v2 lalapps_makeInvetofile lalapps_extractSFTband lalapps_FindSh \
	lalapps_monteint lalapps_skynet getMetric getMesh lalapps_zellepolka	\
	testLatticeCovering

EXTRA_PROGRAMS = $(FRAMEPROGRAMS)


lalapps_zellepolka_SOURCES = zellepolka.c
lalapps_monteint_SOURCES = monteint.c

CFSsources = 	ComputeFStatistic.h	\
		clusters.c clusters.h  	\
		DopplerScan.c DopplerScan.h NeighboringGridPoints.c

lalapps_ComputeFStatistic_SOURCES = ComputeFStatistic.c CFSLALDemod.c $(CFSsources)

lalapps_FrComputeFStatistic_SOURCES = FrComputeFStatistic.c FrComputeFStatistic.h	\
	rngmed.c rngmed.h clusters.c clusters.h  			\
	DopplerScan.c DopplerScan.h NeighboringGridPoints.c

lalapps_FstatShapeTest_SOURCES = FstatShapeTest.c 
lalapps_FstatShapeTestLAL_SOURCES = FstatShapeTestLAL.c

lalapps_SemiAnalyticF_SOURCES =	SemiAnalyticF.c

lalapps_polka_SOURCES = polka.c
lalapps_uberpolka_SOURCES = uberpolka.c


lalapps_makeInvetofile_SOURCES=makeInvetofile.c

lalapps_extractSFTband_SOURCES=extractSFTband.c

lalapps_FindSh_SOURCES=FindSh.c FindSh.h clusters.c clusters.h

lalapps_skynet_SOURCES = skynet.c

compareFstats_SOURCES = compareFstats.c

ComputeFStatistic_v2_SOURCES = ComputeFStatistic_v2.c clusters.c clusters.h  \
	DopplerScan.c DopplerScan.h NeighboringGridPoints.c

TESTS = testFStat.sh

getMetric_SOURCES = getMetric.c DopplerScan.c DopplerScan.h

getMesh_SOURCES = getMesh.c DopplerScan.c DopplerScan.h 

testLatticeCovering_SOURCES = LatticeCovering.c LatticeCovering.h testLatticeCovering.c



