## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libinject.la
libinject_la_SOURCES =	\
	GeneratePPNInspiral.c	SimulateCoherentGW.c \
	InjectVector.c		InjectTimeSeries.c \
	CelestialCoordinates.c	TerrestrialCoordinates.c
## Hack to get this working for now... change later
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include -I$(srcdir)/../../date/include
#INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
