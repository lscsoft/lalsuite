## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libvectorops.la
EXTRA_DIST = MatrixDivide.m4      MatrixMultiplyBase.m4  MiscMatlab.m4 \
MatrixDivideBase.m4  MatrixPower.m4         MiscMatlabBase.m4 \
MatrixMultiply.m4    MatrixPowerBase.m4 \
VectorIndexRangeC.m4 VectorIndexRangeBaseC.m4
MATRIXSRCS = MatrixDivide.c MatrixMultiply.c MatrixPower.c MiscMatlab.c 
libvectorops_la_SOURCES = VectorMultiply.c VectorPolar.c $(MATRIXSRCS) VectorIndexRange.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MAINTAINERCLEANFILES = $(MATRIXSRCS) VectorIndexRange.c
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep

MatrixDivide.c: MatrixDivide.m4 MatrixDivideBase.m4
	cd $(srcdir) && @M4@ MatrixDivide.m4 > MatrixDivide.c

MatrixMultiply.c: MatrixMultiply.m4 MatrixMultiplyBase.m4
	cd $(srcdir) && @M4@ MatrixMultiply.m4 > MatrixMultiply.c

MatrixPower.c: MatrixPower.m4 MatrixPowerBase.m4
	cd $(srcdir) && @M4@ MatrixPower.m4 > MatrixPower.c

MiscMatlab.c: MiscMatlab.m4 MiscMatlabBase.m4
	cd $(srcdir) && @M4@ MiscMatlab.m4 > MiscMatlab.c

VectorIndexRange.c: VectorIndexRangeC.m4 VectorIndexRangeBaseC.m4
	cd $(srcdir) && @M4@ VectorIndexRangeC.m4 > VectorIndexRange.c
@DVIMKFILE@
