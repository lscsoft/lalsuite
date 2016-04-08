include $(top_srcdir)/gnuscripts/lalapps.am

bin_PROGRAMS = \
	lalapps_ComputeAntennaPattern \
	lalapps_FstatMetric_v2 \
	lalapps_PrintDetectorState \
	$(END_OF_LIST)

lalapps_ComputeAntennaPattern_SOURCES = ComputeAntennaPattern.c
lalapps_FstatMetric_v2_SOURCES = FstatMetric_v2.c
lalapps_PrintDetectorState_SOURCES = PrintDetectorState.c

EXTRA_DIST = \
	testComputeAntennaPattern.sh \
	$(END_OF_LIST)

TESTS = \
	testComputeAntennaPattern.sh \
	$(END_OF_LIST)
