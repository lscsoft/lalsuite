## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = Inject.h SimulateCoherentGW.h SimulateSB.h GeneratePPNInspiral.h \
	GenerateTaylorCW.h GenerateSpinOrbitCW.h SkyCoordinates.h SimulateInspiral.h \
	SimulatePopcorn.h GenerateBurst.h GenerateInspiral.h StochasticMC.h \
	GeneratePulsarSignal.h SkyCoordinatesLaldoc.h

MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
all-local: include-link
dvi-local: .dvi-dep
@INCLUDEMKFILE@
@DVIMKFILE@
