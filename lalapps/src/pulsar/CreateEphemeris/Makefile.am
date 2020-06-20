BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

bin_PROGRAMS = lalapps_create_solar_system_ephemeris \
	lalapps_create_time_correction_ephemeris


lalapps_create_solar_system_ephemeris_SOURCES = create_solar_system_ephemeris.c
lalapps_create_time_correction_ephemeris_SOURCES = create_time_correction_ephemeris.c \
	create_time_correction_ephemeris.h

if HAVE_PYTHON
pybin_scripts = lalapps_create_solar_system_ephemeris_python
endif
