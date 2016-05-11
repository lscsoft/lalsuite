include $(top_srcdir)/gnuscripts/lalapps.am

TEST_EXTENSIONS = .sh
SH_LOG_COMPILER = /bin/bash $(srcdir)/test-compiler.sh

if CFITSIO
bin_PROGRAMS = \
	lalapps_WeaveSetup \
	lalapps_Weave \
	$(END_OF_LIST)
else
SH_LOG_COMPILER += skip
endif

lalapps_WeaveSetup_SOURCES = \
	SetupData.c \
	Weave.h \
	WeaveSetup.c \
	$(END_OF_LIST)

lalapps_Weave_SOURCES = \
	CacheResults.c \
	ComputeResults.c \
	OutputData.c \
	SetupData.c \
	Weave.c \
	Weave.h \
	$(END_OF_LIST)

TESTS = \
	TestSingleSegment.sh \
	$(END_OF_LIST)
