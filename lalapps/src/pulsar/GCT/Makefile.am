## Process this file with automake to produce Makefile.in
# $Id$

CFSsrc = $(top_srcdir)/src/pulsar/FDS_isolated
OPTsrc = $(top_srcdir)/src/pulsar/FDS_isolated/OptimizedCFS
EAHsrc = $(top_srcdir)/src/pulsar/hough/src2/EinsteinAtHome

# for using lalapps-defines + functions
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(CFSsrc) -I$(EAHsrc) -I$(top_srcdir)/misc
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

## make sure LALAPPSgitID.h is re-generated at every 'make'
.PHONY: gitID
gitID:
	cd $(top_builddir)/misc && $(MAKE) gitID
BUILT_SOURCES = gitID

bin_PROGRAMS = lalapps_HierarchSearchGCT
EXTRA_PROGRAMS = eah_HierarchSearchGCT eah_HierarchSearchGCT_sse eah_HierarchSearchGCT_sse2

EAHSources = $(EAHsrc)/hs_boinc_extras.c $(EAHsrc)/hs_boinc_options.cpp $(OPTsrc)/LocalComputeFstat.c
HSGCTSources = HierarchSearchGCT.c HierarchSearchGCT.h \
                          GCTtoplist.c GCTtoplist.h \
                          $(CFSsrc)/HeapToplist.c $(CFSsrc)/HeapToplist.h \
                          $(CFSsrc)/ComputeFStatistic.h $(CFSsrc)/Fstat_v3.c \
			  ComputeFstat_RS.c ComputeFstat_RS.h

lalapps_HierarchSearchGCT_SOURCES = $(HSGCTSources)

eah_HierarchSearchGCT_SOURCES = $(HSGCTSources) $(EAHSources)
eah_HierarchSearchGCT_CPPFLAGS = $(AM_CPPFLAGS) -DEAH_BOINC -DHS_OPTIMIZATION

eah_HierarchSearchGCT_sse_SOURCES = $(HSGCTSources) $(EAHSources)
eah_HierarchSearchGCT_sse_CPPFLAGS = $(AM_CPPFLAGS) -DEAH_BOINC -DHS_OPTIMIZATION
eah_HierarchSearchGCT_sse_CFLAGS = $(AM_CFLAGS) -msse -mfpmath=sse

eah_HierarchSearchGCT_sse2_SOURCES = $(HSGCTSources) $(EAHSources)
eah_HierarchSearchGCT_sse2_CPPFLAGS = $(AM_CPPFLAGS) -DEAH_BOINC -DHS_OPTIMIZATION
eah_HierarchSearchGCT_sse2_CFLAGS = $(AM_CFLAGS) -msse -msse2 -mfpmath=sse
