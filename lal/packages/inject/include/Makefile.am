## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = Inject.h SimulateCoherentGW.h GeneratePPNInspiral.h SkyCoordinates.h

all-local: $(top_builddir)/include/lal/Inject.h \
           $(top_builddir)/include/lal/SimulateCoherentGW.h \
           $(top_builddir)/include/lal/GeneratePPNInspiral.h \
           $(top_builddir)/include/lal/SkyCoordinates.h

$(top_builddir)/include/lal/Inject.h: Inject.h
	cp -f $(srcdir)/Inject.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/SimulateCoherentGW.h: SimulateCoherentGW.h
	cp -f $(srcdir)/SimulateCoherentGW.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/GeneratePPNInspiral.h: GeneratePPNInspiral.h
	cp -f $(srcdir)/GeneratePPNInspiral.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/SkyCoordinates.h: SkyCoordinates.h
	cp -f $(srcdir)/SkyCoordinates.h $(top_builddir)/include/lal
