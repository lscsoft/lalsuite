BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalapps_pulsar_test.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

bin_PROGRAMS = \
	lalapps_WeaveSetup \
	lalapps_Weave \
	lalapps_WeaveCompare \
	$(END_OF_LIST)

lalapps_WeaveSetup_SOURCES = \
	SetupData.c \
	SetupData.h \
	Weave.h \
	WeaveSetup.c \
	$(END_OF_LIST)

lalapps_Weave_SOURCES = \
	CacheResults.c \
	CacheResults.h \
	ComputeResults.c \
	ComputeResults.h \
	OutputResults.c \
	OutputResults.h \
	ResultsToplist.c \
	ResultsToplist.h \
	SearchIteration.c \
	SearchIteration.h \
	SearchTiming.c \
	SearchTiming.h \
	SetupData.c \
	SetupData.h \
	Statistics.c \
	Statistics.h \
	Weave.c \
	Weave.h \
	$(END_OF_LIST)

lalapps_WeaveCompare_SOURCES = \
	CacheResults.c \
	OutputResults.c \
	OutputResults.h \
	ResultsToplist.c \
	ResultsToplist.h \
	SearchTiming.c \
	SearchTiming.h \
	SetupData.c \
	SetupData.h \
	Statistics.c \
	Statistics.h \
	ComputeResults.c \
	ComputeResults.h \
	Weave.h \
	WeaveCompare.c \
	$(END_OF_LIST)

# Add shell test scripts to this variable
test_scripts += TestInterpolating.sh
test_scripts += TestNonInterpolating.sh
test_scripts += TestSingleSegment.sh
test_scripts += TestSimulateSearch.sh
test_scripts += TestCacheMaxSize.sh
test_scripts += TestCheckpointing.sh
test_scripts += TestPartitioning.sh

# Add any helper programs required by tests to this variable
test_helpers +=

# All tests require CFITSIO
if !CFITSIO
skip_tests += $(test_scripts)
endif
