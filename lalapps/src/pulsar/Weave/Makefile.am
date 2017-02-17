include $(top_srcdir)/gnuscripts/lalapps.am

TEST_EXTENSIONS = .sh
SH_LOG_COMPILER = /bin/bash $(srcdir)/test-compiler.sh
if !CFITSIO
SH_LOG_COMPILER += skip
endif

bin_PROGRAMS = \
	lalapps_WeaveSetup \
	lalapps_Weave \
	lalapps_WeaveCompare \
	$(END_OF_LIST)

EXTRA_PROGRAMS = \
	lalapps_Weave_Callgrind \
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
	SetupData.c \
	SetupData.h \
	Weave.c \
	Weave.h \
	$(END_OF_LIST)

lalapps_WeaveCompare_SOURCES = \
	OutputResults.c \
	OutputResults.h \
	ResultsToplist.c \
	ResultsToplist.h \
	SetupData.c \
	SetupData.h \
	Weave.h \
	WeaveCompare.c \
	$(END_OF_LIST)

lalapps_Weave_Callgrind_SOURCES = $(lalapps_Weave_SOURCES)
lalapps_Weave_Callgrind_CPPFLAGS = $(AM_CPPFLAGS) -DMAKEFILE_WEAVE_CALLGRIND

TESTS = \
	TestInterpolating.sh \
	TestNonInterpolating.sh \
	TestSingleSegment.sh \
	TestSimulateSearch.sh \
	TestCacheMaxSize.sh \
	TestCheckpointing.sh \
	TestFreqPartition.sh \
	$(END_OF_LIST)

EXTRA_DIST = \
	$(TESTS) \
	test-compiler.sh \
	timestamps-1.txt \
	timestamps-2.txt \
	$(END_OF_LIST)

clean-local:
	-rm -rf Test*.testdir/
