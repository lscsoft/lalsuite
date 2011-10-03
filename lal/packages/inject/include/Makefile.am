pkginclude_HEADERS = \
	Inject.h \
	SimulateCoherentGW.h \
	SkyCoordinates.h

CLEANFILES =
MOSTLYCLEANFILES =
DISTCLEANFILES =

all-local: include-link

@INCLUDEMKFILE@

swig_headers = $(pkginclude_HEADERS)
CLEANFILES += $(swig_clean)
@SWIG_HEADER_MK@
