## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = ComplexFFT.h RealFFT.h

all-local: \
  $(top_builddir)/include/ComplexFFT.h \
  $(top_builddir)/include/RealFFT.h

$(top_builddir)/include/ComplexFFT.h: ComplexFFT.h
	cp -f $(srcdir)/ComplexFFT.h $(top_builddir)/include

$(top_builddir)/include/RealFFT.h: RealFFT.h
	cp -f $(srcdir)/RealFFT.h $(top_builddir)/include
