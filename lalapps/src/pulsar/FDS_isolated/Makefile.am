## Process this file with automake to produce Makefile.in
# $Id$


AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

SUBDIRS = OptimizedCFS ResamplingFstat

bin_PROGRAMS = lalapps_ComputeFStatistic lalapps_FstatShapeTest lalapps_SemiAnalyticF \
	lalapps_polka lalapps_uberpolka lalapps_FstatShapeTestLAL lalapps_compareFstats\
	lalapps_makeInvetofile lalapps_extractSFTband lalapps_FindSh \
	lalapps_monteint lalapps_zellepolka		\
	lalapps_PredictFStat  lalapps_ComputeFStatistic_v2 \
	lalapps_ComputeFStatAnalyticMonteCarloUpperLimit	\
	lalapps_synthesizeBstatMC lalapps_synthesizeTransientStats

check_PROGRAMS = testFstat_v3

EXTRA_PROGRAMS = lalapps_CFSv3 eah_PredictFStat \
		 lalapps_ComputeFStatistic_v2_sse lalapps_ComputeFStatistic_v2_sse2 lalapps_ComputeFStatistic_v2_cuda


# CFSv2 doesn't compile with Condor
#if CONDOR_ENABLED
#EXTRA_PROGRAMS += lalapps_ComputeFStatistic_v2
#else
#bin_PROGRAMS += lalapps_ComputeFStatistic_v2
#endif

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

ComputeFStatisticV2Sources = ComputeFStatistic_v2.c  HeapToplist.c HeapToplist.h OptimizedCFS/ComputeFstatREAL4.c OptimizedCFS/ComputeFstatREAL4.h \
				../transientCW_utils.c ../transientCW_utils.h ../ProbabilityDensity.c ../ProbabilityDensity.h
lalapps_ComputeFStatistic_v2_SOURCES = $(ComputeFStatisticV2Sources)

## optimized ComputeFStatistic_v2 targets

# compilation of SSE targets requires a valid CPU architecture being set in AM_CFLAGS
# (e.g. configure with CFLAGS=-march=pentium-m),
lalapps_ComputeFStatistic_v2_sse_SOURCES = $(ComputeFStatisticV2Sources)
lalapps_ComputeFStatistic_v2_sse_CFLAGS = $(AM_CFLAGS) -msse -mfpmath=sse
lalapps_ComputeFStatistic_v2_sse2_SOURCES = $(ComputeFStatisticV2Sources)
lalapps_ComputeFStatistic_v2_sse2_CFLAGS = $(AM_CFLAGS) -msse -msse2 -mfpmath=sse

# compilation of CUDA target requires cuda header and library directories included in
# CPPFLAGS and LDFLAGS respectively (until LAL/LALAapps have a proper CUDA detection in configure)
lalapps_ComputeFStatistic_v2_cuda_SOURCES = $(ComputeFStatisticV2Sources) OptimizedCFS/ComputeFstatREAL4CUDA.c OptimizedCFS/FStatCUDA.cu
lalapps_ComputeFStatistic_v2_cuda_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_CUDA
lalapps_ComputeFStatistic_v2_cuda_LDADD = $(LDADD) -lcudart

lalapps_PredictFStat_SOURCES = PredictFStat.c ../transientCW_utils.c ../transientCW_utils.h ../ProbabilityDensity.c ../ProbabilityDensity.h
# force linking in C++ mode for E@H
eah_PredictFStat_SOURCES = $(lalapps_PredictFStat_SOURCES) ../hough/src2/EinsteinAtHome/win_lib.cpp

lalapps_ComputeFStatAnalyticMonteCarloUpperLimit_SOURCES = ComputeFStatAnalyticMonteCarloUpperLimit.c

lalapps_synthesizeBstatMC_SOURCES = synthesizeBstatMC.c
lalapps_synthesizeTransientStats_SOURCES = synthesizeTransientStats.c ../transientCW_utils.c ../transientCW_utils.h ../ProbabilityDensity.c ../ProbabilityDensity.h

lalapps_CFSv3_SOURCES = CFSv3.c Fstat_v3.c Fstat_v3.h

testFstat_v3_SOURCES = testFstat_v3.c Fstat_v3.c Fstat_v3.h


EXTRA_DIST = testFStat.sh testCFSv2.sh testGridv2.sh Fstats.ref1 Fstats.ref2 testPredictFStat.sh testCFSv2GPU.sh


TESTS = testFStat.sh testPredictFStat.sh testFstat_v3

if !CONDOR_ENABLED
TESTS += testCFSv2.sh testGridv2.sh ## testCFSv2GPU.sh
endif


TESTS_ENVIRONMENT = \
LALSUITE_TOP_SRCDIR=$(LALSUITE_TOP_SRCDIR) \
LAL_DATA_PATH=.:$${LALPULSAR_PREFIX:+"$${LALPULSAR_PREFIX}/share/lalpulsar"}:$${LALSUITE_TOP_SRCDIR:+"$${LALSUITE_TOP_SRCDIR}/lalpulsar/test/"}

CLEANFILES = testFstat_v3-LFT.sft testFstat_v3-timeseries.dat

# how to compile CUDA files (.cu)
NVCC ?= nvcc
.cu.o:
	$(NVCC) $(CPPFLAGS) -c $< -o $@

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/misc/vcsID.common
