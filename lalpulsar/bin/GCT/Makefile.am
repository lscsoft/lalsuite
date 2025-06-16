BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
man1_MANS =
include $(top_srcdir)/gnuscripts/lalpulsar_bin.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

bin_PROGRAMS = \
	lalpulsar_HierarchSearchGCT \
	$(END_OF_LIST)

EXTRA_PROGRAMS = \
	lalpulsar_HierarchSearchGCT_no_num_count \
	$(END_OF_LIST)

lalpulsar_HierarchSearchGCT_SOURCES = \
	GCTtoplist.c \
	GCTtoplist.h \
	HierarchSearchGCT.c \
	HierarchSearchGCT.h \
	RecalcToplistStats.c \
	RecalcToplistStats.h \
	$(END_OF_LIST)

lalpulsar_HierarchSearchGCT_CPPFLAGS = $(AM_CPPFLAGS)
lalpulsar_HierarchSearchGCT_CFLAGS = $(AM_CFLAGS)

lalpulsar_HierarchSearchGCT_no_num_count_SOURCES = $(lalpulsar_HierarchSearchGCT_SOURCES)
lalpulsar_HierarchSearchGCT_no_num_count_CPPFLAGS = $(AM_CPPFLAGS) -DEXP_NO_NUM_COUNT
lalpulsar_HierarchSearchGCT_no_num_count_CFLAGS = $(AM_CFLAGS)

if HAVE_SSE2_COMPILER

lalpulsar_HierarchSearchGCT_SOURCES += \
	gc_hotloop_sse2.h \
	$(END_OF_LIST)

lalpulsar_HierarchSearchGCT_CPPFLAGS += -DHS_OPTIMIZATION -DHIERARCHSEARCHGCT -DGC_SSE2_OPT
lalpulsar_HierarchSearchGCT_CFLAGS += $(SSE2_FLAGS)

lalpulsar_HierarchSearchGCT_no_num_count_CPPFLAGS += -DHS_OPTIMIZATION -DHIERARCHSEARCHGCT -DGC_SSE2_OPT
lalpulsar_HierarchSearchGCT_no_num_count_CFLAGS += $(SSE2_FLAGS)

endif

# Add shell test scripts to this variable
test_scripts += testHierarchSearchGCT.sh
test_scripts += testHierarchSearchGCT_inject.sh

# Add any helper programs required by tests to this variable
test_helpers +=

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=
