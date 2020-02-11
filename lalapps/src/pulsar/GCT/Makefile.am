BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalapps_pulsar_test.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

AM_CPPFLAGS += -I$(top_srcdir)/src/pulsar/Fstatistic -I$(top_srcdir)/src/pulsar/HoughFstat

bin_PROGRAMS = lalapps_HierarchSearchGCT

EXTRA_PROGRAMS = lalapps_HierarchSearchGCT_SSE2 lalapps_HierarchSearchGCT_SSE2_NONC

lalapps_HierarchSearchGCT_SOURCES = \
	GCTtoplist.c \
	GCTtoplist.h \
	HeapToplist.c \
	HierarchSearchGCT.c \
	HierarchSearchGCT.h \
	RecalcToplistStats.c \
	RecalcToplistStats.h \
	$(END_OF_LIST)

lalapps_HierarchSearchGCT_SSE2_SOURCES = $(HSGCTSources)
lalapps_HierarchSearchGCT_SSE2_CPPFLAGS = $(AM_CPPFLAGS) -DHS_OPTIMIZATION -DHIERARCHSEARCHGCT -DGC_SSE2_OPT
lalapps_HierarchSearchGCT_SSE2_CFLAGS = $(AM_CFLAGS) -msse -msse2 -mfpmath=sse

lalapps_HierarchSearchGCT_SSE2_NONC_SOURCES = $(HSGCTSources)
lalapps_HierarchSearchGCT_SSE2_NONC_CPPFLAGS = $(AM_CPPFLAGS) -DHS_OPTIMIZATION -DHIERARCHSEARCHGCT -DGC_SSE2_OPT -DEXP_NO_NUM_COUNT
lalapps_HierarchSearchGCT_SSE2_NONC_CFLAGS = $(AM_CFLAGS) -msse -msse2 -mfpmath=sse

# Add shell test scripts to this variable
test_scripts += testHierarchSearchGCT.sh
test_scripts += testHierarchSearchGCT_inject.sh

# Add any helper programs required by tests to this variable
test_helpers +=
