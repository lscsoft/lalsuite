## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libfft.la
if INTELFFT
libfft_la_SOURCES = IntelComplexFFT.c IntelRealFFT.c TimeFreqFFT.c
EXTRA_TEXSOURCES = ComplexFFT.c RealFFT.c FFTWMutex.c
else
libfft_la_SOURCES = ComplexFFT.c RealFFT.c FFTWMutex.c TimeFreqFFT.c
EXTRA_TEXSOURCES = 
endif
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
