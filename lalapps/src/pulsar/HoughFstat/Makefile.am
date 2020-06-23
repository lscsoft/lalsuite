BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalapps_pulsar_test.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

AM_CPPFLAGS += -I$(top_srcdir)/src/pulsar/Fstatistic

bin_PROGRAMS = lalapps_HierarchicalSearch
EXTRA_PROGRAMS = lalapps_HierarchicalSearch_sse lalapps_HierarchicalSearch_sse2

lalapps_HierarchicalSearch_SOURCES = \
	HeapToplist.c \
	HierarchicalSearch.c \
	HierarchicalSearch.h \
	HoughFstatToplist.c \
	HoughFstatToplist.h \
	StackSlideFstat.c \
	StackSlideFstat.h \
	$(END_OF_LIST)

# compilation of SSE targets requires a valid CPU architecture being set in AM_CFLAGS
# (e.g. configure with CFLAGS=-march=pentium-m),

lalapps_HierarchicalSearch_sse_SOURCES = $(HierarchicalSearchSources) LocalComputeFstatHoughMap.c
lalapps_HierarchicalSearch_sse_CPPFLAGS = $(AM_CPPFLAGS) -DCOMPUTEFSTATHOUGHMAP=LocalComputeFstatHoughMap
lalapps_HierarchicalSearch_sse_CFLAGS = $(AM_CFLAGS) -msse -mfpmath=sse

lalapps_HierarchicalSearch_sse2_SOURCES = $(HierarchicalSearchSources) LocalComputeFstatHoughMap.c
lalapps_HierarchicalSearch_sse2_CPPFLAGS = $(AM_CPPFLAGS) -DCOMPUTEFSTATHOUGHMAP=LocalComputeFstatHoughMap
lalapps_HierarchicalSearch_sse2_CFLAGS = $(AM_CFLAGS) -msse -msse2 -mfpmath=sse
