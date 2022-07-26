BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps_pulsar.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

bin_PROGRAMS = \
	lalapps_create_solar_system_ephemeris \
	lalapps_create_time_correction_ephemeris \
	$(END_OF_LIST)

lalapps_create_solar_system_ephemeris_SOURCES = \
	create_solar_system_ephemeris.c \
	$(END_OF_LIST)

lalapps_create_time_correction_ephemeris_SOURCES = \
	create_time_correction_ephemeris.c \
	create_time_correction_ephemeris.h \
	$(END_OF_LIST)

if HAVE_PYTHON

pybin_scripts = \
	lalapps_create_solar_system_ephemeris_python \
	$(END_OF_LIST)

endif

# Add shell test scripts to this variable
test_scripts += test_create_solar_system_ephemeris.sh
test_scripts += test_create_solar_system_ephemeris_python.sh
test_scripts += test_create_time_correction_ephemeris.sh

# Add any helper programs required by tests to this variable
test_helpers +=

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=
