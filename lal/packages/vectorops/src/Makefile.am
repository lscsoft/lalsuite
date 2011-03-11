## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libvectorops.la
EXTRA_DIST = \
	MiscMatlab.m4 \
	MiscMatlabBase.m4 \
	VectorIndexRangeC.m4 \
	VectorIndexRangeBaseC.m4

libvectorops_la_SOURCES = \
	MatrixDivide.c \
	MatrixMultiply.c \
	MatrixPower.c \
	MiscMatlab.c \
	VectorMultiply.c \
	VectorPolar.c \
	VectorIndexRange.c

noinst_HEADERS = \
	MatrixDivide_source.c \
	MatrixMultiply_source.c \
	MatrixPower_source.c

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

MAINTAINERCLEANFILES = \
	MiscMatlab.c \
	VectorIndexRange.c

MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex

dvi-local: .dvi-dep

MiscMatlab.c: MiscMatlab.m4 MiscMatlabBase.m4
	cd $(srcdir) && @M4@ MiscMatlab.m4 > MiscMatlab.c

VectorIndexRange.c: VectorIndexRangeC.m4 VectorIndexRangeBaseC.m4
	cd $(srcdir) && @M4@ VectorIndexRangeC.m4 > VectorIndexRange.c
@DVIMKFILE@
