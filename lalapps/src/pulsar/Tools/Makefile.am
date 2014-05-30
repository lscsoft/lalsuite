## Process this file with automake to produce Makefile.in

LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = \
	lalapps_ComputePSD \
	lalapps_GetSFTVelocities \
	lalapps_PrintDetectorState \
	lalapps_ComputeAntennaPattern

lalapps_ComputePSD_SOURCES = ComputePSD.c
lalapps_GetSFTVelocities_SOURCES = GetSFTVelocities.c
lalapps_PrintDetectorState_SOURCES = PrintDetectorState.c
lalapps_ComputeAntennaPattern_SOURCES = ComputeAntennaPattern.c

EXTRA_DIST = \
	testComputePSD.sh \
	testComputeAntennaPattern.sh

TESTS = \
	testComputePSD.sh \
	testComputeAntennaPattern.sh

TESTS_ENVIRONMENT = \
	LAL_DEBUG_LEVEL=msglvl1,memdbg; export LAL_DEBUG_LEVEL; \
	LAL_DATA_PATH="$(LAL_DATA_PATH)"; export LAL_DATA_PATH;

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common
