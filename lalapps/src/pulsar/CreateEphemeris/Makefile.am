## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(top_builddir)/src/lalapps -DPKG_DATA_DIR='"$(pkgdatadir)/"'
LDADD = $(top_builddir)/src/lalapps/liblalapps.la


bin_PROGRAMS = lalapps_create_solar_system_ephemeris \
	       lalapps_create_time_correction_ephemeris


lalapps_create_solar_system_ephemeris_SOURCES = create_solar_system_ephemeris.c
lalapps_create_time_correction_ephemeris_SOURCES = create_time_correction_ephemeris.c create_time_correction_ephemeris.h
