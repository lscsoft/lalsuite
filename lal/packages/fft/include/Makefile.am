## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = ComplexFFT.h RealFFT.h FFTWMutex.h TimeFreqFFT.h
if CUDA
pkginclude_HEADERS += CudaPlan.h
endif
CLEANFILES =
MOSTLYCLEANFILES =
DISTCLEANFILES =
all-local: include-link
@INCLUDEMKFILE@

swig_headers = $(pkginclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
