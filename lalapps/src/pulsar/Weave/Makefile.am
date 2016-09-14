include $(top_srcdir)/gnuscripts/lalapps.am

TEST_EXTENSIONS = .sh
SH_LOG_COMPILER = /bin/bash $(srcdir)/test-compiler.sh

if CFITSIO
bin_PROGRAMS = \
	lalapps_WeaveSetup \
	lalapps_Weave \
	lalapps_WeaveCompare \
	$(END_OF_LIST)
else
SH_LOG_COMPILER += skip
endif

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
	SetupData.c \
	SetupData.h \
	Weave.c \
	Weave.h \
	$(END_OF_LIST)

lalapps_WeaveCompare_SOURCES = \
	OutputResults.c \
	OutputResults.h \
	SetupData.c \
	SetupData.h \
	Weave.h \
	WeaveCompare.c \
	$(END_OF_LIST)

TESTS = \
	TestSingleSegment.sh \
	TestNonInterpolating.sh \
	TestInterpolating.sh \
	TestFreqPartition.sh \
	TestCacheMaxSize.sh \
	TestCheckpointing.sh \
	$(END_OF_LIST)

EXTRA_DIST = \
	$(TESTS) \
	test-compiler.sh \
	timestamps-irregular.txt \
	timestamps-regular.txt \
	$(END_OF_LIST)

clean-local:
	-rm -rf Test*.testdir/
