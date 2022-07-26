BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps_pulsar.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

bin_PROGRAMS = \
	lalapps_HierarchicalSearch \
	$(END_OF_LIST)

lalapps_HierarchicalSearch_SOURCES = \
	HierarchicalSearch.c \
	HierarchicalSearch.h \
	HoughFstatToplist.c \
	HoughFstatToplist.h \
	StackSlideFstat.c \
	StackSlideFstat.h \
	$(END_OF_LIST)

lalapps_HierarchicalSearch_CPPFLAGS = $(AM_CPPFLAGS)
lalapps_HierarchicalSearch_CFLAGS = $(AM_CFLAGS)

if HAVE_SSE2_COMPILER

lalapps_HierarchicalSearch_SOURCES += \
	LocalComputeFstatHoughMap.c \
	hough_sse2.i \
	hough_x64.i \
	hough_x87.i \
	$(END_OF_LIST)

lalapps_HierarchicalSearch_CPPFLAGS += -DCOMPUTEFSTATHOUGHMAP=LocalComputeFstatHoughMap
lalapps_HierarchicalSearch_CFLAGS += $(SSE2_FLAGS)

endif

# Add shell test scripts to this variable
if HAVE_PYTHON
test_scripts += testHierarchicalSearch.py
endif

# Add any helper programs required by tests to this variable
test_helpers +=

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=
