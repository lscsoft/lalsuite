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
test_scripts += testWeave_interpolating.sh
test_scripts += testWeave_non_interpolating.sh
test_scripts += testWeave_single_segment.sh
test_scripts += testWeave_reference_results.sh
test_scripts += testWeave_simulate_search.sh
test_scripts += testWeave_cache_max_size.sh
test_scripts += testWeave_checkpointing.sh
test_scripts += testWeave_partitioning.sh

# Add any helper programs required by tests to this variable
test_helpers +=

# All tests require CFITSIO
if !CFITSIO
skip_tests += $(test_scripts)
endif

# testWeave_reference_results.sh requires output from tests that compare against reference results
testWeave_reference_results.log: testWeave_interpolating.log testWeave_non_interpolating.log testWeave_single_segment.log
