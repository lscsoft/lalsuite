## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = LALStdlib.h LALConstants.h\
        LALAtomicDatatypes.h LALDatatypes.h \
        LALRCSID.h LALStatusMacros.h

EXTRA_DIST = LALPrimer.h

all-local: $(top_builddir)/include/lal/LALStdlib.h \
           $(top_builddir)/include/lal/LALConstants.h \
           $(top_builddir)/include/lal/LALAtomicDatatypes.h \
           $(top_builddir)/include/lal/LALDatatypes.h \
           $(top_builddir)/include/lal/LALRCSID.h \
           $(top_builddir)/include/lal/LALStatusMacros.h

$(top_builddir)/include/lal/LALStdlib.h: LALStdlib.h
	cp -f $(srcdir)/LALStdlib.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/LALConstants.h: LALConstants.h
	cp -f $(srcdir)/LALConstants.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/LALAtomicDatatypes.h: LALAtomicDatatypes.h
	cp -f $(srcdir)/LALAtomicDatatypes.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/LALDatatypes.h: LALDatatypes.h
	cp -f $(srcdir)/LALDatatypes.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/LALRCSID.h: LALRCSID.h
	cp -f $(srcdir)/LALRCSID.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/LALStatusMacros.h: LALStatusMacros.h
	cp -f $(srcdir)/LALStatusMacros.h $(top_builddir)/include/lal

