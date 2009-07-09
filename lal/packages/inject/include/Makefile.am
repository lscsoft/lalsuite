## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = Inject.h SimulateCoherentGW.h GeneratePPNInspiral.h \
	GenerateTaylorCW.h GenerateSpinOrbitCW.h SkyCoordinates.h SimulateInspiral.h \
	GenerateBurst.h GenerateInspiral.h \
	GeneratePulsarSignal.h SkyCoordinatesLaldoc.h NRWaveInject.h \
	InspiralInjectionParams.h GenerateInspRing.h \
	LALSimulation.h LALSimBurst.h LALSimInspiral.h

MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
all-local: include-link
dvi-local: .dvi-dep
@INCLUDEMKFILE@
@DVIMKFILE@
