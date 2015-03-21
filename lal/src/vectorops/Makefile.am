include $(top_srcdir)/gnuscripts/lalsuite_header_links.am

pkginclude_HEADERS = \
	VectorMath.h \
	VectorOps.h \
	$(END_OF_LIST)

noinst_LTLIBRARIES = libvectorops.la

libvectorops_la_SOURCES = \
	VectorMath.c \
	VectorMath_AVX.c \
	VectorMath_FPU.c \
	VectorMath_SSE.c \
	VectorMath_avx_mathfun.h \
	VectorMath_sse_mathfun.h \
	VectorMultiply.c \
	VectorPolar.c \
	$(END_OF_LIST)

noinst_HEADERS = \
	$(END_OF_LIST)
