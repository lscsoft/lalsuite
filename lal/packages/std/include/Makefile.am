## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = LALStdlib.h LALConstants.h\
        LALAtomicDatatypes.h LALDatatypes.h \
        LALRCSID.h LALStatusMacros.h \
	LALStdio.h LALVersion.h LALMalloc.h LALErrno.h LALError.h LALGSL.h \
	XLALError.h XLALGSL.h \
	StringInput.h Grid.h
EXTRA_DIST = LALPrimer.h
EXTRA_TEXSOURCES = LALPrimer.h
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
all-local: include-link
dvi-local: .dvi-dep
@INCLUDEMKFILE@
@DVIMKFILE@
