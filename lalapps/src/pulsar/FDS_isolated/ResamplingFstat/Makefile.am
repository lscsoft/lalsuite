## Process this file with automake to produce Makefile.in
# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

## make sure VCSInfo.h is re-generated at every 'make'
.PHONY: gitID
gitID:
	cd $(top_builddir)/src/lalapps && $(MAKE) liblalapps.la
BUILT_SOURCES = gitID

bin_PROGRAMS = lalapps_ComputeFStatistic_resamp

EXTRA_PROGRAMS = resampling

resampling_SOURCES = proc.c fstat.c

lalapps_ComputeFStatistic_resamp_SOURCES = ComputeFStatistic_resamp.c ../HeapToplist.c ../HeapToplist.h

TESTS = testCFS_resamp.sh
EXTRA_DIST = testCFS_resamp.sh


TESTS_ENVIRONMENT = \
LAL_TOP_SRCDIR=$(LAL_TOP_SRCDIR) \
LAL_DATA_PATH=.:$${LALPULSAR_PREFIX:+"$${LALPULSAR_PREFIX}/share/lalpulsar"}:$${LAL_TOP_SRCDIR:+"$${LAL_TOP_SRCDIR}/../lalpulsar/test/"}
