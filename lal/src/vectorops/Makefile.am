include $(top_srcdir)/gnuscripts/lalsuite_header_links.am

pkginclude_HEADERS = \
	VectorMath.h \
	VectorOps.h \
	$(END_OF_LIST)

noinst_LTLIBRARIES = libvectorops.la

libvectorops_la_SOURCES = \
	VectorMath.c \
	VectorMath_GEN.c \
	VectorMath_GEN_Find.c \
	VectorMultiply.c \
	VectorPolar.c \
	$(END_OF_LIST)

noinst_HEADERS = \
	VectorMath_avx_mathfun.h \
	VectorMath_internal.h \
	VectorMath_sse_mathfun.h \
	$(END_OF_LIST)

libvectorops_la_LIBADD =

if HAVE_SSE_COMPILER
noinst_LTLIBRARIES += libvectormath_sse.la
libvectorops_la_LIBADD += libvectormath_sse.la
libvectormath_sse_la_SOURCES = VectorMath_SSEx.c
libvectormath_sse_la_CFLAGS = $(AM_CFLAGS) $(SSE_CFLAGS)
endif

if HAVE_SSE2_COMPILER
noinst_LTLIBRARIES += libvectormath_sse2.la
libvectorops_la_LIBADD += libvectormath_sse2.la
libvectormath_sse2_la_SOURCES = VectorMath_SSEx.c
libvectormath_sse2_la_CFLAGS = $(AM_CFLAGS) $(SSE2_CFLAGS)
endif

if HAVE_SSSE3_COMPILER
noinst_LTLIBRARIES += libvectormath_ssse3.la
libvectorops_la_LIBADD += libvectormath_ssse3.la
libvectormath_ssse3_la_SOURCES = VectorMath_SSSE3_Find.c
libvectormath_ssse3_la_CFLAGS = $(AM_CFLAGS) $(SSSE3_CFLAGS)
endif

if HAVE_AVX_COMPILER
noinst_LTLIBRARIES += libvectormath_avx.la
libvectorops_la_LIBADD += libvectormath_avx.la
libvectormath_avx_la_SOURCES = VectorMath_AVXx.c
libvectormath_avx_la_CFLAGS = $(AM_CFLAGS) $(AVX_CFLAGS)
endif

if HAVE_AVX2_COMPILER
noinst_LTLIBRARIES += libvectormath_avx2.la
libvectorops_la_LIBADD += libvectormath_avx2.la
libvectormath_avx2_la_SOURCES = VectorMath_AVXx.c VectorMath_AVX2_Find.c
libvectormath_avx2_la_CFLAGS = $(AM_CFLAGS) $(AVX2_CFLAGS)
endif
