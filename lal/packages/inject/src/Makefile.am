noinst_LTLIBRARIES = libinject.la
libinject_la_SOURCES = \
	CelestialCoordinates.c \
	CelestialCoordinatesLaldoc.c \
	InjectTimeSeries.c \
	InjectVector.c \
	LALSimulation.c \
	SimulateCoherentGW.c \
	SkyCoordinates.c \
	TerrestrialCoordinates.c \
	check_series_macros.h

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex

dvi-local: .dvi-dep

@DVIMKFILE@
