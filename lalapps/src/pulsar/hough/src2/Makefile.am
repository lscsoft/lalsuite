 ## Process this file with automake to produce Makefile.in
# $Id$

# for using lalapps-defines + functions
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(top_srcdir)/src/pulsar/FDS_isolated -I$(top_srcdir)/misc

## Compiling with OpenCL on Linux
# LDFLAGS = -lOpenCL
#
## Compiling with OpenCL on Mac
# LDFLAGS = -framework OpenCL

bin_PROGRAMS = lalapps_HoughValidateAM lalapps_HoughMismatch lalapps_HierarchicalSearch

EXTRA_PROGRAMS = lalapps_HierarchicalSearch_v2 lalapps_HierarchicalSearchGC \
		 lalapps_HierarchicalSearch_sse lalapps_HierarchicalSearch_sse2 lalapps_HierarchicalSearch_cuda \
		 eah_HierarchicalSearch eah_HierarchicalSearch_sse eah_HierarchicalSearch_sse2 eah_HierarchicalSearch_cuda		 

## make sure LALAPPSgitID.h is re-generated at every 'make'
.PHONY: gitID
gitID:
	cd $(top_builddir)/misc && $(MAKE) gitID
BUILT_SOURCES = gitID

## make sure lalapps.la is built
lalappsla = $(top_builddir)/src/lalapps/liblalapps.la
$(lalappsla):
	cd $(top_builddir)/src/lalapps && $(MAKE)
LDADD = $(lalappsla)

# how to compile CUDA files (.cu)
NVCC ?= nvcc
.cu.o:
	$(NVCC) $(CPPFLAGS) -c $< -o $@

FDSsrc = ../../FDS_isolated

lalapps_HoughValidateAM_SOURCES = HoughValidateAM.c ../src/MCInjectHoughS2.h \
	                  ../src/SFTbin.h ../src/SFTbin.c \
			  ../src/DriveHoughColor.h ../src/PeakSelect.h \
			  ../src/PeakSelect.c

lalapps_HoughMismatch_SOURCES = HoughMismatch.c ../src/MCInjectHoughS2.h \
			../src/SFTbin.h ../src/SFTbin.c \
			../src/DriveHoughColor.h ../src/PeakSelect.h \
			../src/PeakSelect.c

lalapps_HierarchicalSearch_v2_SOURCES = HierarchicalSearch_v2.c HierarchicalSearch.h StackSlideFstat.c StackSlideFstat.h \
			  $(FDSsrc)/FstatToplist.c $(FDSsrc)/FstatToplist.h \
			  $(FDSsrc)/HeapToplist.c $(FDSsrc)/HeapToplist.h \
			  $(FDSsrc)/ComputeFStatistic.h

lalapps_HierarchicalSearchGC_SOURCES = HierarchicalSearchGC.c HierarchicalSearchGC.h \
			  StackSlideFstat.c StackSlideFstat.h \
			  HoughFStatToplist.c HoughFStatToplist.h \
			  $(FDSsrc)/HeapToplist.c $(FDSsrc)/HeapToplist.h \
			  $(FDSsrc)/ComputeFStatistic.h

HierarchicalSearchSources = HierarchicalSearch.c HierarchicalSearch.h \
			  StackSlideFstat.c StackSlideFstat.h \
			  HoughFStatToplist.c HoughFStatToplist.h \
			  $(FDSsrc)/HeapToplist.c $(FDSsrc)/HeapToplist.h \
			  $(FDSsrc)/ComputeFStatistic.h \
			  $(FDSsrc)/OptimizedCFS/ComputeFstatREAL4.h $(FDSsrc)/OptimizedCFS/ComputeFstatREAL4.c

lalapps_HierarchicalSearch_SOURCES = $(HierarchicalSearchSources)


## optimized HierarchicalSearch targets

# compilation of SSE targets requires a valid CPU architecture being set in AM_CFLAGS
# (e.g. configure with CFLAGS=-march=pentium-m),
lalapps_HierarchicalSearch_sse_SOURCES = $(HierarchicalSearchSources)
lalapps_HierarchicalSearch_sse_CFLAGS = $(AM_CFLAGS) -msse -mfpmath=sse
lalapps_HierarchicalSearch_sse2_SOURCES = $(HierarchicalSearchSources)
lalapps_HierarchicalSearch_sse2_CFLAGS = $(AM_CFLAGS) -msse -msse2 -mfpmath=sse

# compilation of CUDA target requires cuda header and library directories included in
# CPPFLAGS and LDFLAGS respectively (until LAL/LALAapps have a proper CUDA detection in configure)
lalapps_HierarchicalSearch_cuda_SOURCES = $(HierarchicalSearchSources) $(FDSsrc)/OptimizedCFS/ComputeFstatREAL4CUDA.c $(FDSsrc)/OptimizedCFS/FStatCUDA.cu
lalapps_HierarchicalSearch_cuda_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_CUDA
lalapps_HierarchicalSearch_cuda_LDADD = $(LDADD) -lcudart

# Einstein@home targets
EAHSources = $(FDSsrc)/OptimizedCFS/LocalComputeFstat.c EinsteinAtHome/LocalComputeFstatHoughMap.c EinsteinAtHome/hs_boinc_extras.c EinsteinAtHome/hs_boinc_options.cpp

eah_HierarchicalSearch_SOURCES = $(HierarchicalSearchSources) $(EAHSources)
eah_HierarchicalSearch_CPPFLAGS = $(AM_CPPFLAGS) -DEAH_BOINC -DHS_OPTIMIZATION

eah_HierarchicalSearch_sse_SOURCES = $(HierarchicalSearchSources) $(EAHSources)
eah_HierarchicalSearch_sse_CPPFLAGS = $(AM_CPPFLAGS) -DEAH_BOINC -DHS_OPTIMIZATION
eah_HierarchicalSearch_sse_CFLAGS = $(AM_CFLAGS) -msse -mfpmath=sse

eah_HierarchicalSearch_sse2_SOURCES = $(HierarchicalSearchSources) $(EAHSources)
eah_HierarchicalSearch_sse2_CPPFLAGS = $(AM_CPPFLAGS) -DEAH_BOINC -DHS_OPTIMIZATION
eah_HierarchicalSearch_sse2_CFLAGS = $(AM_CFLAGS) -msse -msse2 -mfpmath=sse

eah_HierarchicalSearch_cuda_SOURCES = $(HierarchicalSearchSources) $(EAHSources) $(FDSsrc)/OptimizedCFS/ComputeFstatREAL4CUDA.c $(FDSsrc)/OptimizedCFS/FStatCUDA.cu
eah_HierarchicalSearch_cuda_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_CUDA -DEAH_BOINC -DHS_OPTIMIZATION
