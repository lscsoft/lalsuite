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
	XLALGSL.h

if SWIG_BUILD
pkginclude_HEADERS += \
	swig-lal-error.i
endif

EXTRA_DIST = LALPrimer.h

EXTRA_TEXSOURCES = LALPrimer.h

CLEANFILES =
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex

all-local: include-link

dvi-local: .dvi-dep

@INCLUDEMKFILE@
@DVIMKFILE@

swig_headers = $(pkginclude_HEADERS)
include $(top_srcdir)/swig/swig-header.am
