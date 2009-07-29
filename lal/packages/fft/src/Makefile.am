## Process this file with automake to produce Makefile.in
## $Id$
if QTHREAD
QTHREADSRC = qthread.c
else
QTHREADSRC =
endif

if INTELFFT
FFTSRC = IntelComplexFFT.c IntelRealFFT.c $(QTHREADSRC)
EXTRA_TEXSOURCES = ComplexFFT.c RealFFT.c FFTWMutex.c
else
FFTSRC = ComplexFFT.c RealFFT.c FFTWMutex.c
EXTRA_TEXSOURCES =
endif

noinst_LTLIBRARIES = libfft.la
libfft_la_SOURCES = $(FFTSRC) TimeFreqFFT.c AverageSpectrum.c Convolution.c

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
EXTRA_DIST = qthread.c IntelComplexFFT.c IntelRealFFT.c \
             ComplexFFT.c RealFFT.c FFTWMutex.c TimeFreqFFT.c
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
