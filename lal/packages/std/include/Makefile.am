## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = LALConfig.h LALStdlib.h LALConstants.h\
        LALAtomicDatatypes.h LALDatatypes.h \
        LALRCSID.h LALStatusMacros.h

EXTRA_DIST = LALPrimer.h

all-local: $(top_builddir)/include/LALConfig.h \
           $(top_builddir)/include/LALStdlib.h \
           $(top_builddir)/include/LALConstants.h \
           $(top_builddir)/include/LALAtomicDatatypes.h \
           $(top_builddir)/include/LALDatatypes.h \
           $(top_builddir)/include/LALRCSID.h \
           $(top_builddir)/include/LALStatusMacros.h

$(top_builddir)/include/LALConfig.h: LALConfig.h
	cp -f $(srcdir)/LALConfig.h $(top_builddir)/include

$(top_builddir)/include/LALStdlib.h: LALStdlib.h
	cp -f $(srcdir)/LALStdlib.h $(top_builddir)/include

$(top_builddir)/include/LALConstants.h: LALConstants.h
	cp -f $(srcdir)/LALConstants.h $(top_builddir)/include

$(top_builddir)/include/LALAtomicDatatypes.h: LALAtomicDatatypes.h
	cp -f $(srcdir)/LALAtomicDatatypes.h $(top_builddir)/include

$(top_builddir)/include/LALDatatypes.h: LALDatatypes.h
	cp -f $(srcdir)/LALDatatypes.h $(top_builddir)/include

$(top_builddir)/include/LALRCSID.h: LALRCSID.h
	cp -f $(srcdir)/LALRCSID.h $(top_builddir)/include

$(top_builddir)/include/LALStatusMacros.h: LALStatusMacros.h
	cp -f $(srcdir)/LALStatusMacros.h $(top_builddir)/include

