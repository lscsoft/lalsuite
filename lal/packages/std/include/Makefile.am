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
	LALRCSID.h \
	LALStatusMacros.h \
	LALStdio.h \
	LALStdlib.h \
	LALString.h \
	LALVersion.h \
	StringInput.h \
	XLALError.h \
	XLALGSL.h

EXTRA_DIST = LALPrimer.h

EXTRA_TEXSOURCES = LALPrimer.h

MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex

all-local: include-link

dvi-local: .dvi-dep

@INCLUDEMKFILE@
@DVIMKFILE@
