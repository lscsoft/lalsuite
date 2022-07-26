BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps_pulsar.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

bin_PROGRAMS = \
	lalapps_HierarchSearchGCT \
	$(END_OF_LIST)

EXTRA_PROGRAMS = \
	lalapps_HierarchSearchGCT_NONC \
	$(END_OF_LIST)

lalapps_HierarchSearchGCT_SOURCES = \
	GCTtoplist.c \
	GCTtoplist.h \
	HierarchSearchGCT.c \
	HierarchSearchGCT.h \
	RecalcToplistStats.c \
	RecalcToplistStats.h \
	$(END_OF_LIST)

lalapps_HierarchSearchGCT_CPPFLAGS = $(AM_CPPFLAGS)
lalapps_HierarchSearchGCT_CFLAGS = $(AM_CFLAGS)

lalapps_HierarchSearchGCT_NONC_SOURCES = $(lalapps_HierarchSearchGCT_SOURCES)
lalapps_HierarchSearchGCT_NONC_CPPFLAGS = $(AM_CPPFLAGS) -DEXP_NO_NUM_COUNT
lalapps_HierarchSearchGCT_NONC_CFLAGS = $(AM_CFLAGS)

if HAVE_SSE2_COMPILER

lalapps_HierarchSearchGCT_SOURCES += \
	gc_hotloop_sse2.h \
	$(END_OF_LIST)

lalapps_HierarchSearchGCT_CPPFLAGS += -DHS_OPTIMIZATION -DHIERARCHSEARCHGCT -DGC_SSE2_OPT
lalapps_HierarchSearchGCT_CFLAGS += $(SSE2_FLAGS)

lalapps_HierarchSearchGCT_NONC_CPPFLAGS += -DHS_OPTIMIZATION -DHIERARCHSEARCHGCT -DGC_SSE2_OPT
lalapps_HierarchSearchGCT_NONC_CFLAGS += $(SSE2_FLAGS)

endif

# Add shell test scripts to this variable
test_scripts += testHierarchSearchGCT.sh
test_scripts += testHierarchSearchGCT_inject.sh

# Add any helper programs required by tests to this variable
test_helpers +=

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=
