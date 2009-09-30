## Process this file with automake to produce Makefile.in
# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(top_srcdir)/misc
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

SUBDIRS = ResamplingFstat

## make sure LALAPPSgitID.h is re-generated at every 'make'
.PHONY: gitID
gitID:
	cd $(top_builddir)/misc && $(MAKE) gitID
BUILT_SOURCES = gitID

bin_PROGRAMS = lalapps_ComputeFStatistic lalapps_FstatShapeTest lalapps_SemiAnalyticF \
	lalapps_polka lalapps_uberpolka lalapps_FstatShapeTestLAL lalapps_compareFstats\
	lalapps_makeInvetofile lalapps_extractSFTband lalapps_FindSh \
	lalapps_monteint lalapps_zellepolka		\
	lalapps_ComputeFStatistic_v2 lalapps_PredictFStat \
	lalapps_ComputeFStatAnalyticMonteCarloUpperLimit	\
	lalapps_synthesizeBstatMC

lalapps_zellepolka_SOURCES = zellepolka.c
lalapps_monteint_SOURCES = monteint.c

CFSsources = 	ComputeFStatistic.h	\
		clusters.c clusters.h  	\
		FstatToplist.c FstatToplist.h \
		HeapToplist.c HeapToplist.h

lalapps_ComputeFStatistic_SOURCES = ComputeFStatistic.c CFSLALDemod.c $(CFSsources)

##lalapps_FrComputeFStatistic_SOURCES = FrComputeFStatistic.c FrComputeFStatistic.h	\
##	rngmed.c rngmed.h clusters.c clusters.h  			\
##	DopplerScan.c DopplerScan.h NeighboringGridPoints.c

lalapps_FstatShapeTest_SOURCES = FstatShapeTest.c
lalapps_FstatShapeTestLAL_SOURCES = FstatShapeTestLAL.c

lalapps_SemiAnalyticF_SOURCES =	SemiAnalyticF.c

lalapps_polka_SOURCES = polka.c
lalapps_uberpolka_SOURCES = uberpolka.c

lalapps_makeInvetofile_SOURCES=makeInvetofile.c

lalapps_extractSFTband_SOURCES=extractSFTband.c

lalapps_FindSh_SOURCES=FindSh.c FindSh.h clusters.c clusters.h

lalapps_compareFstats_SOURCES = compareFstats.c

lalapps_ComputeFStatistic_v2_SOURCES = 	ComputeFStatistic_v2.c  HeapToplist.c HeapToplist.h

lalapps_PredictFStat_SOURCES = PredictFStat.c

lalapps_ComputeFStatAnalyticMonteCarloUpperLimit_SOURCES = ComputeFStatAnalyticMonteCarloUpperLimit.c

lalapps_synthesizeBstatMC_SOURCES = synthesizeBstatMC.c

EXTRA_DIST = testFStat.sh testCFSv2.sh testGridv2.sh Fstats.ref1 Fstats.ref2 testPredictFStat.sh

TESTS = testFStat.sh testCFSv2.sh testGridv2.sh testPredictFStat.sh
TESTS_ENVIRONMENT = \
LAL_TOP_SRCDIR=$(LAL_TOP_SRCDIR) \
LAL_PREFIX=$(LAL_PREFIX) \
LAL_DATA_PATH=.:$${LAL_TOP_SRCDIR:+"$${LAL_TOP_SRCDIR}/packages/pulsar/test"}:$${LAL_PREFIX:+"$${LAL_PREFIX}/share/lal"}

CLEANFILES =
