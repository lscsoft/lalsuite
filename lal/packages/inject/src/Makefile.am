## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libinject.la
libinject_la_SOURCES =	\
	GetInspiralParams.c		GeneratePPNInspiral.c \
	GeneratePPNAmpCorInspiral.c \
	GenerateSpinOrbitCW.c		GenerateEllipticSpinOrbitCW.c \
	GenerateParabolicSpinOrbitCW.c	GenerateHyperbolicSpinOrbitCW.c \
	GenerateTaylorCW.c		SimulateCoherentGW.c \
	InjectVector.c			InjectTimeSeries.c \
	CelestialCoordinates.c		TerrestrialCoordinates.c \
	SkyCoordinates.c		SimulateInspiral.c \
	GenerateInspiral.c		GenerateInspiralSmooth.c  \
	GeneratePulsarSignal.c		SimulatePulsarSignal.c	 \
	CelestialCoordinatesLaldoc.c    NRWaveInject.c \
	InspiralInjectionParams.c    GenerateInspRing.c	 \
	LALSimulation.c \
	LALSimInspiral.c 		LALSimInspiralPNMode.c \
	check_series_macros.h

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
