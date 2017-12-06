BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/lalapps.am

bin_PROGRAMS =  lalapps_skynet \
		lalapps_getMetric \
		lalapps_getMesh \
		lalapps_FstatMetric_v2

lalapps_skynet_SOURCES = skynet.c

lalapps_getMetric_SOURCES = getMetric.c

lalapps_getMesh_SOURCES = getMesh.c

lalapps_FstatMetric_v2_SOURCES = FstatMetric_v2.c


EXTRA_DIST =

TESTS =

TESTS_ENVIRONMENT = \
	LAL_DEBUG_LEVEL=msglvl1,memdbg; export LAL_DEBUG_LEVEL; \
	LAL_DATA_PATH="$(LAL_DATA_PATH)"; export LAL_DATA_PATH;

CLEANFILES =
