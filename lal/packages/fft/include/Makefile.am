## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = ComplexFFT.h RealFFT.h

all-local: \
  $(top_srcdir)/include/ComplexFFT.h \
  $(top_srcdir)/include/RealFFT.h

$(top_srcdir)/include/ComplexFFT.h: ComplexFFT.h
	cp -f ComplexFFT.h $(top_srcdir)/include

$(top_srcdir)/include/RealFFT.h: RealFFT.h
	cp -f RealFFT.h $(top_srcdir)/include
