include $(top_srcdir)/gnuscripts/lalsuite_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_cuda.am

if CUDA
CUDAHDRS = CudaPlan.h
endif

pkginclude_HEADERS = \
	ComplexFFT.h \
	RealFFT.h \
	FFTWMutex.h \
	TimeFreqFFT.h \
	$(CUDAHDRS)

if QTHREAD
QTHREADSRC = qthread.c
else
QTHREADSRC =
endif

if INTELFFT
FFTSRC = \
	IntelComplexFFT.c \
	IntelRealFFT.c \
	FFTWMutex.c \
	$(QTHREADSRC)
FFTHDR =
FFTCXXSRC =
FFTCXXGENSRC =
FFTLIBCXX =
else
if CUDA
FFTSRC = \
	CudaComplexFFT.c \
	CudaRealFFT.c \
	FFTWMutex.c \
	CudaFunctions.c \
	$(END_OF_LIST)
FFTHDR =
FFTCXXSRC =
FFTCXXGENSRC = CudaFFT.cpp
FFTLIBCXX = libfftcxx.la
else
FFTSRC = \
	ComplexFFT.c \
	RealFFT.c \
	FFTWMutex.c \
	$(END_OF_LIST)
FFTHDR = \
	RealFFT_source.c \
	ComplexFFT_source.c \
	$(END_OF_LIST)
FFTCXXSRC =
FFTCXXGENSRC =
FFTLIBCXX =
endif
endif

noinst_LTLIBRARIES = \
	libfft.la \
	$(FFTLIBCXX)

libfft_la_SOURCES = \
	TimeFreqFFT.c \
	AverageSpectrum.c \
	Convolution.c \
	$(FFTSRC)

noinst_HEADERS = \
	$(FFTHDR)

libfft_la_LIBADD = $(FFTLIBCXX)
libfftcxx_la_SOURCES = $(FFTCXXSRC)
nodist_libfftcxx_la_SOURCES = $(FFTCXXGENSRC)

CLEANFILES = CudaFFT.cpp

EXTRA_DIST = \
	ComplexFFT.c \
	CudaComplexFFT.c \
	CudaFFT.cu \
	CudaFFT.h \
	CudaFunctions.h \
	CudaRealFFT.c \
	FFTWMutex.c \
	IntelComplexFFT.c \
	IntelRealFFT.c \
	RealFFT.c \
	RealFFT_source.c \
	TimeFreqFFT.c \
	qthread.c \
	$(END_OF_LIST)
