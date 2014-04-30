## Process this file with automake to produce Makefile.in


AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -DPKG_DATA_DIR='"$(pkgdatadir)/"' @AM_CPPFLAGS@
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

SUBDIRS = ResamplingFstat

bin_PROGRAMS = \
	lalapps_ComputeFStatistic \
	lalapps_FstatShapeTest \
	lalapps_SemiAnalyticF \
	lalapps_polka \
	lalapps_uberpolka \
	lalapps_FstatShapeTestLAL \
	lalapps_compareFstats\
	lalapps_makeInvetofile \
	lalapps_extractSFTband \
	lalapps_FindSh \
	lalapps_zellepolka \
	lalapps_PredictFStat \
	lalapps_ComputeFStatistic_v2 \
	lalapps_ComputeFStatAnalyticMonteCarloUpperLimit \
	lalapps_synthesizeBstatMC \
	lalapps_synthesizeTransientStats \
	lalapps_synthesizeLVStats

EXTRA_PROGRAMS = \
	lalapps_CFSv3 \
	eah_PredictFStat \
	eah_ComputeFStatistic_v2

lalapps_zellepolka_SOURCES = zellepolka.c

CFSsources = \
	ComputeFStatistic.h \
	clusters.c \
	clusters.h \
	FstatToplist.c \
	FstatToplist.h \
	HeapToplist.c \
	HeapToplist.h

lalapps_ComputeFStatistic_SOURCES = \
	ComputeFStatistic.c \
	CFSLALDemod.c \
	$(CFSsources)

lalapps_FstatShapeTest_SOURCES = FstatShapeTest.c

lalapps_FstatShapeTestLAL_SOURCES = FstatShapeTestLAL.c

lalapps_SemiAnalyticF_SOURCES = SemiAnalyticF.c

lalapps_polka_SOURCES = polka.c

lalapps_uberpolka_SOURCES = uberpolka.c

lalapps_makeInvetofile_SOURCES=makeInvetofile.c

lalapps_extractSFTband_SOURCES=extractSFTband.c

lalapps_FindSh_SOURCES=FindSh.c FindSh.h clusters.c clusters.h

lalapps_compareFstats_SOURCES = compareFstats.c

ComputeFStatisticV2Sources = \
	ComputeFStatistic_v2.c \
	HeapToplist.c \
	HeapToplist.h

lalapps_ComputeFStatistic_v2_SOURCES = $(ComputeFStatisticV2Sources)

lalapps_PredictFStat_SOURCES = PredictFStat.c

# force linking in C++ mode for E@H
eah_PredictFStat_SOURCES = \
	$(lalapps_PredictFStat_SOURCES) \
	../hough/src2/EinsteinAtHome/win_lib.cpp
eah_ComputeFStatistic_v2_SOURCES = \
	$(ComputeFStatisticV2Sources) \
	../hough/src2/EinsteinAtHome/win_lib.cpp

lalapps_ComputeFStatAnalyticMonteCarloUpperLimit_SOURCES = ComputeFStatAnalyticMonteCarloUpperLimit.c

lalapps_synthesizeBstatMC_SOURCES = synthesizeBstatMC.c

lalapps_synthesizeTransientStats_SOURCES = synthesizeTransientStats.c

lalapps_synthesizeLVStats_SOURCES = synthesizeLVStats.c

lalapps_CFSv3_SOURCES = CFSv3.c

EXTRA_DIST = \
	Fstats.ref1 \
	Fstats.ref2 \
	testPredictFStat.sh \
	testCFSv2.sh \
	testCFSv2_grids.sh \
	testCFSv2_resamp.sh \
	testCFSv2_binary.sh

EXTRA_DIST += \
	testCFSv2.dat.ref.gz \
	testCFSv2_grid0.dat.ref.gz \
	testCFSv2_grid1.dat.ref.gz \
	testCFSv2_grid2.dat.ref.gz \
	testCFSv2_grid6.dat.ref.gz \
	testCFSv2_grid8.dat.ref.gz \
	testCFSv2_grid9.dat.ref.gz

TESTS = testPredictFStat.sh

if !CONDOR_ENABLED
TESTS += \
	testCFSv2.sh \
	testCFSv2_grids.sh \
	testCFSv2_resamp.sh \
	testCFSv2_binary.sh
endif

TESTS_ENVIRONMENT = \
	LAL_DEBUG_LEVEL=msglvl1,memdbg; export LAL_DEBUG_LEVEL; \
	LAL_DATA_PATH="$(LAL_DATA_PATH)"; export LAL_DATA_PATH;

# how to compile CUDA files (.cu)
NVCC ?= nvcc
.cu.o:
	$(NVCC) $(CPPFLAGS) -c $< -o $@

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common
