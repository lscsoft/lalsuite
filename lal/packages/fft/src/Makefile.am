## Process this file with automake to produce Makefile.in
if QTHREAD
QTHREADSRC = qthread.c
else
QTHREADSRC =
endif

if INTELFFT
FFTSRC = IntelComplexFFT.c IntelRealFFT.c $(QTHREADSRC)
FFTCXXSRC =
FFTCXXGENSRC =
FFTLIBCXX =
else
if CUDA
FFTSRC = CudaComplexFFT.c CudaRealFFT.c FFTWMutex.c FFTWMutex.h CudaFunctions.c
FFTCXXSRC =
FFTCXXGENSRC = CudaFFT.cpp
FFTLIBCXX = libfftcxx.la
else
FFTSRC = ComplexFFT.c RealFFT.c FFTWMutex.c FFTWMutex.h
FFTCXXSRC =
FFTCXXGENSRC =
FFTLIBCXX =
endif
endif

noinst_LTLIBRARIES = libfft.la $(FFTLIBCXX)
libfft_la_SOURCES = $(FFTSRC) TimeFreqFFT.c AverageSpectrum.c Convolution.c
libfft_la_LIBADD = $(FFTLIBCXX)
libfftcxx_la_SOURCES = $(FFTCXXSRC)
nodist_libfftcxx_la_SOURCES = $(FFTCXXGENSRC)

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
EXTRA_DIST = qthread.c IntelComplexFFT.c IntelRealFFT.c \
             ComplexFFT.c RealFFT.c FFTWMutex.c FFTWMutex.h TimeFreqFFT.c \
             CudaComplexFFT.c CudaRealFFT.c CudaFFT.cu CudaFunctions.h \
             CudaFFT.h

CLEANFILES = CudaFFT.cpp
MOSTLYCLEANFILES =
DISTCLEANFILES =

# cuda build rules
include $(top_srcdir)/misc/cuda.common
