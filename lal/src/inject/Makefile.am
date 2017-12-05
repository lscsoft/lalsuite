include $(top_srcdir)/gnuscripts/lalsuite_header_links.am

pkginclude_HEADERS = \
	SimulateCoherentGW.h \
	SkyCoordinates.h \
	$(END_OF_LIST)

noinst_LTLIBRARIES = libinject.la

libinject_la_SOURCES = \
	CelestialCoordinates.c \
	SimulateCoherentGW.c \
	SkyCoordinates.c \
	TerrestrialCoordinates.c \
	$(END_OF_LIST)

EXTRA_DIST = \
	SkyCoordinates.dox \
	SkyCoordinates_conversions.png \
	inject_ecliptic.png \
	inject_galactic.png \
	inject_geodetic.png \
	inject_geodeticsing.png \
	inject_geographic.png \
	inject_horizon.png \
	inject_lat_long.png \
	inject_phase_diagram.png \
	$(END_OF_LIST)
