if CUDA
CUDAHDRS = CudaPlan.h
endif

pkginclude_HEADERS = \
	ComplexFFT.h \
	RealFFT.h \
	FFTWMutex.h \
	TimeFreqFFT.h \
	$(CUDAHDRS)

all-local: include-link
@INCLUDEMKFILE@

swig_headers = $(pkginclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
