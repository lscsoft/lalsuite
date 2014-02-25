## Process this file with automake to produce Makefile.in

# for using lalapps-defines + functions
AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(top_srcdir)/src/pulsar/FDS_isolated -DPKG_DATA_DIR='"$(pkgdatadir)/"' @AM_CPPFLAGS@
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

## Compiling with OpenCL on Linux
# LDFLAGS = -lOpenCL
#
## Compiling with OpenCL on Mac
# LDFLAGS = -framework OpenCL

# how to compile CUDA files (.cu)
NVCC ?= nvcc
.cu.o:
	$(NVCC) $(CPPFLAGS) -c $< -o $@
# allow to specify a different CUDA runtime library (e.g. for Windows when cross-compiling)
CUDART ?= -lcudart

bin_PROGRAMS = lalapps_HoughValidateAM lalapps_HoughMismatch lalapps_HierarchicalSearch lalapps_RecalcHSCandidates

EXTRA_PROGRAMS = lalapps_HierarchicalSearch_sse lalapps_HierarchicalSearch_sse2 \
		 eah_HierarchicalSearch

# lalapps_HierarchicalSearch doesn't condor compile
#if CONDOR_ENABLED
#EXTRA_PROGRAMS += lalapps_HierarchicalSearch
#else
#bin_PROGRAMS += lalapps_HierarchicalSearch
#endif

FDSsrc = ../../FDS_isolated
Recalcsrc = ../../GCT

lalapps_HoughValidateAM_SOURCES = HoughValidateAM.c ../src/MCInjectHoughS2.h \
	                  ../src/SFTbin.h ../src/SFTbin.c \
			  ../src/DriveHoughColor.h ../src/PeakSelect.h \
			  ../src/PeakSelect.c

lalapps_HoughMismatch_SOURCES = HoughMismatch.c ../src/MCInjectHoughS2.h \
			../src/SFTbin.h ../src/SFTbin.c \
			../src/DriveHoughColor.h ../src/PeakSelect.h \
			../src/PeakSelect.c

HierarchicalSearchSources = HierarchicalSearch.c HierarchicalSearch.h \
			  StackSlideFstat.c StackSlideFstat.h \
			  HoughFStatToplist.c HoughFStatToplist.h \
			  $(FDSsrc)/HeapToplist.c $(FDSsrc)/HeapToplist.h \
			  $(FDSsrc)/ComputeFStatistic.h

lalapps_HierarchicalSearch_SOURCES = $(HierarchicalSearchSources)

lalapps_RecalcHSCandidates_SOURCES = RecalcHSCandidates.c HierarchicalSearch.h \
				  HoughFStatToplist.c HoughFStatToplist.h \
				  $(FDSsrc)/HeapToplist.c $(FDSsrc)/HeapToplist.h \
				  $(FDSsrc)/ComputeFStatistic.h \
				  $(Recalcsrc)/RecalcToplistStats.c $(Recalcsrc)/RecalcToplistStats.h


HSLocalSources = EinsteinAtHome/LocalComputeFstatHoughMap.c
EAHSources = $(HSLocalSources) EinsteinAtHome/hs_boinc_extras.c EinsteinAtHome/hs_boinc_options.cpp EinsteinAtHome/win_lib.cpp

## optimized HierarchicalSearch targets

# compilation of SSE targets requires a valid CPU architecture being set in AM_CFLAGS
# (e.g. configure with CFLAGS=-march=pentium-m),
lalapps_HierarchicalSearch_sse_SOURCES = $(HierarchicalSearchSources) $(HSLocalSources)
lalapps_HierarchicalSearch_sse_CPPFLAGS = $(AM_CPPFLAGS) -DCOMPUTEFSTATHOUGHMAP=LocalComputeFstatHoughMap
lalapps_HierarchicalSearch_sse_CFLAGS = $(AM_CFLAGS) -msse -mfpmath=sse
lalapps_HierarchicalSearch_sse2_SOURCES = $(HierarchicalSearchSources) $(HSLocalSources)
lalapps_HierarchicalSearch_sse2_CPPFLAGS = $(AM_CPPFLAGS) -DCOMPUTEFSTATHOUGHMAP=LocalComputeFstatHoughMap
lalapps_HierarchicalSearch_sse2_CFLAGS = $(AM_CFLAGS) -msse -msse2 -mfpmath=sse

# Einstein@home targets
eah_HierarchicalSearch_SOURCES = $(HierarchicalSearchSources) $(EAHSources)
eah_HierarchicalSearch_CPPFLAGS = $(AM_CPPFLAGS) -DEAH_BOINC -DHS_OPTIMIZATION
# this shouldn't harm - it just prefers libraries in the build directory over e.g. system ones
eah_HierarchicalSearch_LDFLAGS = -L. $(AM_LDFLAGS)
# RELEASE_DEPS contains additional dependencies of the E@H build
# items in there should have rules listed below to create them
# the idea is to increase compatibility of E@H applications by linking these libs statically
eah_HierarchicalSearch_DEPENDENCIES = $(RELEASE_DEPS)
eah_HierarchicalSearch_LDADD = $(LDADD) $(RELEASE_LDADD)

libstdc++.a:
	cp -f `$(CC) -print-file-name=$@` $@

exchndl.o: EinsteinAtHome/exchndl.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

erp_execinfo_plus.o: EinsteinAtHome/erp_execinfo_plus.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

TESTS = testRecalcHS.sh
EXTRA_DIST = testRecalcHS.sh

TESTS_ENVIRONMENT = \
	LAL_DEBUG_LEVEL=msglvl1,memdbg; export LAL_DEBUG_LEVEL; \
	LAL_DATA_PATH="$(LAL_DATA_PATH)"; export LAL_DATA_PATH;

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common
