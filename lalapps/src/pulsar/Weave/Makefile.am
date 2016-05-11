include $(top_srcdir)/gnuscripts/lalapps.am

TEST_EXTENSIONS = .sh
SH_LOG_COMPILER = /bin/bash $(srcdir)/test-compiler.sh

if CFITSIO
bin_PROGRAMS = \
	$(END_OF_LIST)
else
SH_LOG_COMPILER += skip
endif

TESTS = \
	$(END_OF_LIST)
