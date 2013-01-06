EXTRA_DIST = LALConfig.dox \
	LALDatatypes.dox \
	LALError.dox \
	LALMalloc.dox \
	LALStatusMacros.dox

pkginclude_HEADERS = \
	Grid.h \
	LALAtomicDatatypes.h \
	LALComplex.h \
	LALConstants.h \
	LALDatatypes.h \
	LALErrno.h \
	LALError.h \
	LALGSL.h \
	LALMalloc.h \
	LALStatusMacros.h \
	LALStdio.h \
	LALStdlib.h \
	LALString.h \
	LALVersion.h \
	StringInput.h \
	XLALError.h \
	XLALGSL.h \
	lalswig_error.i \
	lalswig_gsl.i \
	lalswig_ligotimegps.i

all-local: include-link
@INCLUDEMKFILE@

swig_headers = $(pkginclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
