## Process this file with automake to produce Makefile.in
## $Id$
noinst_LTLIBRARIES = libfft.la
if QTHREAD
QTHREADSRC = qthread.c
else
QTHREADSRC = 
endif
if INTELFFT
libfft_la_SOURCES = IntelComplexFFT.c IntelRealFFT.c TimeFreqFFT.c AverageSpectrum.c $(QTHREADSRC)
EXTRA_TEXSOURCES = ComplexFFT.c RealFFT.c FFTWMutex.c
else
libfft_la_SOURCES = ComplexFFT.c RealFFT.c FFTWMutex.c TimeFreqFFT.c AverageSpectrum.c
EXTRA_TEXSOURCES = 
endif
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
EXTRA_DIST = qthread.c IntelComplexFFT.c IntelRealFFT.c \
             ComplexFFT.c RealFFT.c FFTWMutex.c TimeFreqFFT.c
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
