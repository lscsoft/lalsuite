include $(top_srcdir)/gnuscripts/lalapps.am

TEST_EXTENSIONS = .sh
SH_LOG_COMPILER = /bin/bash $(srcdir)/test-compiler.sh

if CFITSIO
bin_PROGRAMS = \
	lalapps_WeaveSetup \
	$(END_OF_LIST)
else
SH_LOG_COMPILER += skip
endif

lalapps_WeaveSetup_SOURCES = \
	SetupData.c \
	Weave.h \
	WeaveSetup.c \
	$(END_OF_LIST)

TESTS = \
	$(END_OF_LIST)
