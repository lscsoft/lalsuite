## Process this file with automake to produce Makefile.in
# $Id$
pkginclude_HEADERS = VectorOps.h Matrix.h VectorIndexRange.h
MAINTAINERCLEANFILES = VectorIndexRange.h
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex 
all-local: include-link
dvi-local: .dvi-dep
EXTRA_DIST = VectorIndexRangeBaseH.m4 VectorIndexRangeH.m4

VectorIndexRange.h: VectorIndexRangeBaseH.m4 VectorIndexRangeH.m4
	cd $(srcdir) && @M4@ VectorIndexRangeH.m4 > VectorIndexRange.h
@INCLUDEMKFILE@
@DVIMKFILE@
