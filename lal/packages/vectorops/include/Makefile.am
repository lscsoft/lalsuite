pkginclude_HEADERS = \
	VectorOps.h \
	Matrix.h \
	VectorIndexRange.h

CLEANFILES =
MOSTLYCLEANFILES =
DISTCLEANFILES =
MAINTAINERCLEANFILES =

all-local: include-link

@INCLUDEMKFILE@

swig_headers = $(pkginclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
