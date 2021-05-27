BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalapps_pulsar_test.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

bin_PROGRAMS = \
	lalapps_ComputeAntennaPattern \
	lalapps_FstatMetric_v2 \
	lalapps_PrintDetectorState \
	$(END_OF_LIST)

lalapps_ComputeAntennaPattern_SOURCES = ComputeAntennaPattern.c
lalapps_FstatMetric_v2_SOURCES = FstatMetric_v2.c
lalapps_PrintDetectorState_SOURCES = PrintDetectorState.c

# Add shell test scripts to this variable
test_scripts += testComputeAntennaPattern.sh

# Add any helper programs required by tests to this variable
test_helpers +=

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=
