BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
man1_MANS =
include $(top_srcdir)/gnuscripts/lalpulsar_bin.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

bin_PROGRAMS = \
	lalpulsar_HierarchicalSearch \
	$(END_OF_LIST)

lalpulsar_HierarchicalSearch_SOURCES = \
	HierarchicalSearch.c \
	HierarchicalSearch.h \
	HoughFstatToplist.c \
	HoughFstatToplist.h \
	StackSlideFstat.c \
	StackSlideFstat.h \
	$(END_OF_LIST)

lalpulsar_HierarchicalSearch_CPPFLAGS = $(AM_CPPFLAGS)
lalpulsar_HierarchicalSearch_CFLAGS = $(AM_CFLAGS)

if HAVE_SSE2_COMPILER

lalpulsar_HierarchicalSearch_SOURCES += \
	LocalComputeFstatHoughMap.c \
	hough_sse2.i \
	hough_x64.i \
	hough_x87.i \
	$(END_OF_LIST)

lalpulsar_HierarchicalSearch_CPPFLAGS += -DCOMPUTEFSTATHOUGHMAP=LocalComputeFstatHoughMap
lalpulsar_HierarchicalSearch_CFLAGS += $(SSE2_FLAGS)

endif

# Add shell test scripts to this variable
if HAVE_PYTHON
test_scripts += testHierarchicalSearch.py
endif

# Add any helper programs required by tests to this variable
test_helpers +=

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=
