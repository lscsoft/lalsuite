include $(top_srcdir)/gnuscripts/lalsuite_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_cuda.am

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
	CudaFunctions.c
FFTHDR =
FFTCXXSRC =
FFTCXXGENSRC = CudaFFT.cpp
FFTLIBCXX = libfftcxx.la
else
FFTSRC = \
	ComplexFFT.c \
	RealFFT.c \
	FFTWMutex.c
FFTHDR = \
	RealFFT_source.c \
	ComplexFFT_source.c
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

EXTRA_DIST = \
	qthread.c \
	IntelComplexFFT.c \
	IntelRealFFT.c \
	ComplexFFT.c \
	RealFFT.c \
	RealFFT_source.c \
	FFTWMutex.c \
	TimeFreqFFT.c \
	CudaComplexFFT.c \
	CudaRealFFT.c \
	CudaFFT.cu \
	CudaFunctions.h \
	CudaFFT.h

CLEANFILES = CudaFFT.cpp
