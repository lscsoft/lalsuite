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
if CUDA
FFTSRC = CudaComplexFFT.c CudaRealFFT.c FFTWMutex.c CudaFunctions.c CudaFFT.cu
else
FFTSRC = ComplexFFT.c RealFFT.c FFTWMutex.c
EXTRA_TEXSOURCES =
endif
endif

noinst_LTLIBRARIES = libfft.la
libfft_la_SOURCES = $(FFTSRC) TimeFreqFFT.c AverageSpectrum.c Convolution.c

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
EXTRA_DIST = qthread.c IntelComplexFFT.c IntelRealFFT.c \
             ComplexFFT.c RealFFT.c FFTWMutex.c TimeFreqFFT.c \
             CudaComplexFFT.c CudaRealFFT.c CudaFFT.cu CudaFunctions.h \
             CudaFFT.h
if CUDA
MAINTAINERCLEANFILES = CudaFFT.cu.c CudaFFT.cu.loT
endif
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@

if CUDA
SUFFIXES = .cu .c
.cu.c:
	@nvcc -cuda --host-compilation=c $(INCLUDES) $(CPPFLAGS) --output-file $@ $<
endif
