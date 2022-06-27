BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps_pulsar.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

bin_PROGRAMS = lalapps_HierarchicalSearch
EXTRA_PROGRAMS = lalapps_HierarchicalSearch_sse lalapps_HierarchicalSearch_sse2

lalapps_HierarchicalSearch_SOURCES = \
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

# Add shell test scripts to this variable
if HAVE_PYTHON
test_scripts += testHierarchicalSearch.py
endif
