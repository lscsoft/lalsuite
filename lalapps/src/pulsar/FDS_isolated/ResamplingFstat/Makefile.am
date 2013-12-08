## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = @AM_CPPFLAGS@ -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -DPKG_DATA_DIR='"$(pkgdatadir)/"'
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_ComputeFStatistic_resamp

EXTRA_PROGRAMS = resampling

resampling_SOURCES = proc.c fstat.c

lalapps_ComputeFStatistic_resamp_SOURCES = ComputeFStatistic_resamp.c ../HeapToplist.c ../HeapToplist.h

TESTS = testCFS_resamp.sh
EXTRA_DIST = testCFS_resamp.sh


TESTS_ENVIRONMENT = \
	LAL_DEBUG_LEVEL=memdbg; export LAL_DEBUG_LEVEL; \
	LAL_DATA_PATH="$(LAL_DATA_PATH)"; export LAL_DATA_PATH;

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common
