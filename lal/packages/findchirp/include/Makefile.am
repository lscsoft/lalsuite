## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = FindChirp.h \
                     FindChirpChisq.h \
                     FindChirpSP.h \
                     FindChirpExch.h \
                     FindChirpEngine.h

all-local: \
        $(top_builddir)/include/lal/FindChirp.h \
        $(top_builddir)/include/lal/FindChirpChisq.h \
        $(top_builddir)/include/lal/FindChirpSP.h \
        $(top_builddir)/include/lal/FindChirpExch.h \
        $(top_builddir)/include/lal/FindChirpEngine.h 

$(top_builddir)/include/lal/FindChirp.h: FindChirp.h
	cp -f $(srcdir)/FindChirp.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/FindChirpChisq.h: FindChirpChisq.h
	cp -f $(srcdir)/FindChirpChisq.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/FindChirpSP.h: FindChirpSP.h
	cp -f $(srcdir)/FindChirpSP.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/FindChirpExch.h: FindChirpExch.h
	cp -f $(srcdir)/FindChirpExch.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/FindChirpEngine.h: FindChirpEngine.h
	cp -f $(srcdir)/FindChirpEngine.h $(top_builddir)/include/lal
