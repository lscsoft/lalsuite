## Process this file with automake to produce Makefile.in
# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(top_srcdir)/misc
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_skynet lalapps_getMetric lalapps_testFlatPulsarMetric lalapps_matchGrids \
	lalapps_getMesh lalapps_FstatMetric lalapps_testLatticeCovering  \
	lalapps_FstatMetric_v2

## make sure LALAPPSgitID.h is re-generated at every 'make'
.PHONY: gitID
gitID:
	cd $(top_builddir)/misc && $(MAKE) gitID
BUILT_SOURCES = gitID


lalapps_testLatticeCovering_SOURCES = testLatticeCovering.c

lalapps_testFlatLatticeTilingPulsar_SOURCES = testFlatLatticeTilingPulsar.c

lalapps_skynet_SOURCES = skynet.c

lalapps_getMetric_SOURCES = getMetric.c

lalapps_getMesh_SOURCES = getMesh.c

lalapps_testFlatPulsarMetric_SOURCES = testFlatPulsarMetric.c

lalapps_matchGrids_SOURCES = matchGrids.c

lalapps_FstatMetric_SOURCES = FstatMetric.c

lalapps_FstatMetric_v2_SOURCES = FstatMetric_v2.c
