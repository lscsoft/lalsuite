BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps_pulsar.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

bin_PROGRAMS = \
	lalapps_ComputeAntennaPattern \
	lalapps_FstatMetric_v2 \
	lalapps_PrintDetectorState \
	$(END_OF_LIST)

lalapps_ComputeAntennaPattern_SOURCES = \
	ComputeAntennaPattern.c \
	$(END_OF_LIST)

lalapps_FstatMetric_v2_SOURCES = \
	FstatMetric_v2.c \
	$(END_OF_LIST)

lalapps_PrintDetectorState_SOURCES = \
	PrintDetectorState.c \
	$(END_OF_LIST)

# Add shell test scripts to this variable
test_scripts += testComputeAntennaPattern.sh
test_scripts += testFstatMetric_v2.sh

# Add any helper programs required by tests to this variable
test_helpers +=

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=
