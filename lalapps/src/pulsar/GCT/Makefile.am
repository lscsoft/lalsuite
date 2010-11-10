## Process this file with automake to produce Makefile.in
# $Id$

CFSsrc = $(top_srcdir)/src/pulsar/FDS_isolated
OPTsrc = $(top_srcdir)/src/pulsar/FDS_isolated/OptimizedCFS
EAHsrc = $(top_srcdir)/src/pulsar/hough/src2/EinsteinAtHome

# for using lalapps-defines + functions
INCLUDES = -I$(top_srcdir)/src/lalapps -I$(CFSsrc) -I$(EAHsrc)
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_HierarchSearchGCT
EXTRA_PROGRAMS = eah_HierarchSearchGCT lalapps_HierarchSearchGCT_cuda lalapps_HierarchSearchGCT_VFStat

EAHSources = $(EAHsrc)/hs_boinc_extras.c $(EAHsrc)/hs_boinc_options.cpp $(EAHsrc)/win_lib.cpp $(OPTsrc)/LocalComputeFstat.c
HSGCTSources = HierarchSearchGCT.c HierarchSearchGCT.h \
                          GCTtoplist.c GCTtoplist.h \
                          $(CFSsrc)/HeapToplist.c $(CFSsrc)/HeapToplist.h \
                          $(CFSsrc)/ComputeFStatistic.h $(CFSsrc)/Fstat_v3.c \
			  ComputeFstat_RS.c ComputeFstat_RS.h

lalapps_HierarchSearchGCT_SOURCES = $(HSGCTSources)

lalapps_HierarchSearchGCT_cuda_SOURCES = $(HSGCTSources) $(OPTsrc)/ComputeFstatREAL4.c $(OPTsrc)/ComputeFstatREAL4CUDA.c $(OPTsrc)/FStatCUDA.cu
lalapps_HierarchSearchGCT_cuda_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_CUDA
lalapps_HierarchSearchGCT_cuda_LDADD = $(LDADD) -lcudart

lalapps_HierarchSearchGCT_VFStat_SOURCES = $(HSGCTSources) $(OPTsrc)/LocalComputeFstat.c
lalapps_HierarchSearchGCT_VFStat_CPPFLAGS = $(AM_CPPFLAGS) -DHS_OPTIMIZATION -DHIERARCHSEARCHGCT

eah_HierarchSearchGCT_SOURCES = $(HSGCTSources) $(EAHSources)
eah_HierarchSearchGCT_CPPFLAGS = $(AM_CPPFLAGS) -DEAH_BOINC -DHS_OPTIMIZATION -DHIERARCHSEARCHGCT
# this shouldn't harm - it just prefers libraries in the build directory over e.g. system ones
eah_HierarchSearchGCT_LDFLAGS = -L. $(AM_LDFLAGS)
# RELEASE_DEPS contains additional dependencies of the E@H build
# items in there should have rules listed below to create them
# the idea is to increase compatibility of E@H applications by linking these libs statically
eah_HierarchSearchGCT_DEPENDENCIES = $(RELEASE_DEPS)
eah_HierarchSearchGCT_LDADD = $(LDADD) $(RELEASE_LDADD)

libstdc++.a:
	cp -f `$(CC) -print-file-name=$@` $@
libz.a:
	cp -f `$(CC) -print-file-name=$@` $@

exchndl.o: $(EAHsrc)/exchndl.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

erp_execinfo_plus.o: $(EAHsrc)/erp_execinfo_plus.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/misc/vcsID.common

# how to compile CUDA files (.cu)
NVCC ?= nvcc
.cu.o:
	$(NVCC) $(CPPFLAGS) -c $< -o $@

## pretty dirty hack, particulary for compiling an E@H Mac OS 10.3 App on 10.4
eah_HierarchSearchGCT_manual: eah_HierarchSearchGCT-HierarchSearchGCT.o \
    eah_HierarchSearchGCT-GCTtoplist.o \
    eah_HierarchSearchGCT-HeapToplist.o \
    eah_HierarchSearchGCT-Fstat_v3.o \
    eah_HierarchSearchGCT-ComputeFstat_RS.o \
    eah_HierarchSearchGCT-hs_boinc_extras.o \
    eah_HierarchSearchGCT-hs_boinc_options.o \
    eah_HierarchSearchGCT-win_lib.o \
    eah_HierarchSearchGCT-LocalComputeFstat.o \
    libstdc++.a
	$(CXX) $(CPPFLAGS) $(AC_CPPFLAGS) -o $@ $^ $(LDFLAGS) $(AC_LDFLAGS) $(LDADD) $(AC_LDADD) $(LIBS) $(AC_LIBS) libstdc++.a
