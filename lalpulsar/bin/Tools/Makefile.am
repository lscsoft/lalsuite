BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
man1_MANS =
include $(top_srcdir)/gnuscripts/lalpulsar_bin.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

bin_PROGRAMS = \
	lalpulsar_ComputeAntennaPattern \
	lalpulsar_FstatMetric_v2 \
	lalpulsar_PrintDetectorState \
	$(END_OF_LIST)

lalpulsar_ComputeAntennaPattern_SOURCES = \
	ComputeAntennaPattern.c \
	$(END_OF_LIST)

lalpulsar_FstatMetric_v2_SOURCES = \
	FstatMetric_v2.c \
	$(END_OF_LIST)

lalpulsar_PrintDetectorState_SOURCES = \
	PrintDetectorState.c \
	$(END_OF_LIST)

# Add shell test scripts to this variable
test_scripts += testComputeAntennaPattern.sh
test_scripts += testFstatMetric_v2.sh

# Add any helper programs required by tests to this variable
test_helpers +=

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=
