## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = LALConfig.h LALStdlib.h LALConstants.h\
        LALAtomicDatatypes.h LALDatatypes.h \
        LALRCSID.h LALStatusMacros.h

EXTRA_DIST = LALPrimer.h

all-local: $(top_srcdir)/include/LALConfig.h \
           $(top_srcdir)/include/LALStdlib.h \
           $(top_srcdir)/include/LALConstants.h \
           $(top_srcdir)/include/LALAtomicDatatypes.h \
           $(top_srcdir)/include/LALDatatypes.h \
           $(top_srcdir)/include/LALRCSID.h \
           $(top_srcdir)/include/LALStatusMacros.h

$(top_srcdir)/include/LALConfig.h: LALConfig.h
	cp -f LALConfig.h $(top_srcdir)/include

$(top_srcdir)/include/LALStdlib.h: LALStdlib.h
	cp -f LALStdlib.h $(top_srcdir)/include

$(top_srcdir)/include/LALConstants.h: LALConstants.h
	cp -f LALConstants.h $(top_srcdir)/include

$(top_srcdir)/include/LALAtomicDatatypes.h: LALAtomicDatatypes.h
	cp -f LALAtomicDatatypes.h $(top_srcdir)/include

$(top_srcdir)/include/LALDatatypes.h: LALDatatypes.h
	cp -f LALDatatypes.h $(top_srcdir)/include

$(top_srcdir)/include/LALRCSID.h: LALRCSID.h
	cp -f LALRCSID.h $(top_srcdir)/include

$(top_srcdir)/include/LALStatusMacros.h: LALStatusMacros.h
	cp -f LALStatusMacros.h $(top_srcdir)/include

