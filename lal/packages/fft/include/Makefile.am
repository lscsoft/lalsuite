## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = ComplexFFT.h RealFFT.h

all-local: \
  $(top_builddir)/include/lal/ComplexFFT.h \
  $(top_builddir)/include/lal/RealFFT.h

$(top_builddir)/include/lal/ComplexFFT.h: ComplexFFT.h
	cp -f $(srcdir)/ComplexFFT.h $(top_builddir)/include/lal

$(top_builddir)/include/lal/RealFFT.h: RealFFT.h
	cp -f $(srcdir)/RealFFT.h $(top_builddir)/include/lal
