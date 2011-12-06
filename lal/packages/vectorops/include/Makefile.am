pkginclude_HEADERS = \
	VectorOps.h \
	Matrix.h \
	VectorIndexRange.h

CLEANFILES =
MOSTLYCLEANFILES =
DISTCLEANFILES =
MAINTAINERCLEANFILES = \
	VectorIndexRange.h

all-local: include-link

EXTRA_DIST = \
	VectorIndexRangeBaseH.m4 \
	VectorIndexRangeH.m4 \
	VectorPairTypedefsBase.m4

VectorIndexRange.h: VectorIndexRangeBaseH.m4 VectorIndexRangeH.m4 VectorPairTypedefsBase.m4
	cd $(srcdir) && @M4@ VectorIndexRangeH.m4 > VectorIndexRange.h

@INCLUDEMKFILE@

swig_headers = $(pkginclude_HEADERS)
include $(top_srcdir)/swig/swig-header.am
