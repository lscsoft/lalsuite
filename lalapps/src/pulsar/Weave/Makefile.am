BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps_pulsar.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_cuda.am

bin_PROGRAMS = \
	lalpulsar_Weave \
	lalpulsar_WeaveCompare \
	lalpulsar_WeaveConcat \
	lalpulsar_WeaveSetup \
	$(END_OF_LIST)

lalpulsar_WeaveSetup_SOURCES = \
	SetupData.c \
	SetupData.h \
	Weave.h \
	WeaveSetup.c \
	$(END_OF_LIST)

lalpulsar_Weave_SOURCES = \
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

lalpulsar_WeaveCompare_SOURCES = \
	CacheResults.c \
	CacheResults.h \
	ComputeResults.c \
	ComputeResults.h \
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
	Weave.h \
	WeaveCompare.c \
	$(END_OF_LIST)

lalpulsar_WeaveConcat_SOURCES = \
	CacheResults.c \
	CacheResults.h \
	ComputeResults.c \
	ComputeResults.h \
	OutputResults.c \
	OutputResults.h \
	ResultsToplist.c \
	ResultsToplist.h \
	SearchTiming.c \
	SearchTiming.h \
	Statistics.c \
	Statistics.h \
	Weave.h \
	WeaveConcat.c \
	$(END_OF_LIST)

noinst_LTLIBRARIES =

if CUDA

noinst_LTLIBRARIES += libweave_computeresults_cuda.la
LDADD = libweave_computeresults_cuda.la
libweave_computeresults_cuda_la_SOURCES = ComputeResults_CUDA.cu
libweave_computeresults_cuda_la_LDFLAGS = -no-undefined
MOSTLYCLEANFILES += ComputeResults_CUDA.cpp

endif

# Add shell test scripts to this variable
test_scripts += testWeave_interpolating.sh
test_scripts += testWeave_non_interpolating.sh
test_scripts += testWeave_single_segment.sh
test_scripts += testWeave_reference_results.sh
test_scripts += testWeave_simulate_search.sh
test_scripts += testWeave_cache_max_size.sh
test_scripts += testWeave_checkpointing.sh
test_scripts += testWeave_partitioning.sh
test_scripts += testWeave_histograms.sh
test_scripts += testWeave_concatenation.sh
test_scripts += testWeave_random_injection.sh

# Add any helper programs required by tests to this variable
test_helpers +=

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=

# All tests require CFITSIO
if !CFITSIO

skip_tests += $(test_scripts)

endif

# testWeave_reference_results.sh requires output from tests that compare against reference results
testWeave_reference_results.log: testWeave_interpolating.log testWeave_non_interpolating.log testWeave_single_segment.log
