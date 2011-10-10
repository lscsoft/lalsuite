## Process this file with automake to produce Makefile.in
# $Id$
pkginclude_HEADERS = VectorOps.h Matrix.h VectorIndexRange.h
CLEANFILES =
MAINTAINERCLEANFILES = VectorIndexRange.h
MOSTLYCLEANFILES =
DISTCLEANFILES =
all-local: include-link

EXTRA_DIST = VectorIndexRangeBaseH.m4 VectorIndexRangeH.m4 VectorPairTypedefsBase.m4

VectorIndexRange.h: VectorIndexRangeBaseH.m4 VectorIndexRangeH.m4 VectorPairTypedefsBase.m4
	cd $(srcdir) && @M4@ VectorIndexRangeH.m4 > VectorIndexRange.h
@INCLUDEMKFILE@

swig_headers = $(pkginclude_HEADERS)
CLEANFILES += $(swig_clean)
@SWIG_HEADER_MK@
