BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalpulsar_bin.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

bin_PROGRAMS = \
	lalpulsar_create_time_correction_ephemeris \
	$(END_OF_LIST)

lalpulsar_create_time_correction_ephemeris_SOURCES = \
	create_time_correction_ephemeris.c \
	create_time_correction_ephemeris.h \
	$(END_OF_LIST)

# Add shell test scripts to this variable
test_scripts += test_create_time_correction_ephemeris.sh

# Add any helper programs required by tests to this variable
test_helpers +=

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=
