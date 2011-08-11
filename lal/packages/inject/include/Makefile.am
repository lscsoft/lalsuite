pkginclude_HEADERS = \
	Inject.h \
	SimulateCoherentGW.h \
	SkyCoordinates.h \
	SkyCoordinatesLaldoc.h

CLEANFILES =
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex

all-local: include-link

dvi-local: .dvi-dep

@INCLUDEMKFILE@

@DVIMKFILE@

swig_headers = $(pkginclude_HEADERS)
CLEANFILES += $(swig_clean)
@SWIG_HEADER_MK@
