noinst_LTLIBRARIES = libinject.la
libinject_la_SOURCES = \
	CelestialCoordinates.c \
	CelestialCoordinatesLaldoc.c \
	GenerateEllipticSpinOrbitCW.c \
	GenerateHyperbolicSpinOrbitCW.c \
	GenerateInspRing.c \
	GenerateInspiral.c \
	GenerateInspiralSmooth.c  \
	GeneratePPNAmpCorInspiral.c \
	GeneratePPNInspiral.c \
	GenerateParabolicSpinOrbitCW.c \
	GeneratePulsarSignal.c \
	GenerateSpinOrbitCW.c \
	GenerateTaylorCW.c \
	GetInspiralParams.c \
	InjectTimeSeries.c \
	InjectVector.c \
	InspiralInjectionParams.c \
	LALSimInspiral.c \
	LALSimInspiralPNMode.c \
	LALSimulation.c \
	NRWaveInject.c \
	SimulateCoherentGW.c \
	SimulateInspiral.c \
	SimulatePulsarSignal.c \
	SkyCoordinates.c \
	TerrestrialCoordinates.c \
	check_series_macros.h

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex

dvi-local: .dvi-dep

@DVIMKFILE@
