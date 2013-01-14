## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_skynet lalapps_getMetric lalapps_testFlatPulsarMetric lalapps_matchGrids \
	lalapps_getMesh lalapps_FstatMetric lalapps_testLatticeCovering  \
	lalapps_FstatMetric_v2

lalapps_testLatticeCovering_SOURCES = testLatticeCovering.c

lalapps_skynet_SOURCES = skynet.c

lalapps_getMetric_SOURCES = getMetric.c

lalapps_getMesh_SOURCES = getMesh.c

lalapps_testFlatPulsarMetric_SOURCES = testFlatPulsarMetric.c

lalapps_matchGrids_SOURCES = matchGrids.c

lalapps_FstatMetric_SOURCES = FstatMetric.c

lalapps_FstatMetric_v2_SOURCES = FstatMetric_v2.c


EXTRA_DIST = testMetricCodes.py

TESTS = testMetricCodes.py

TESTS_ENVIRONMENT = LALPULSAR_DATADIR="$(LALPULSAR_DATADIR)"

CLEANFILES = Fmetric.dat Fmetric_v2.dat

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/misc/vcsID.common
