## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libutilities.la
libutilities_la_SOURCES = \
  Random.c MersenneRandom.c FindRoot.c Integrate.c Interpolate.c \
  HeapSort.c ODE.c CoarseGrainFrequencySeries.c Dirichlet.c \
  MatrixOps.c DetInverse.c DetInverseInternal.c Eigen.c EigenInternal.c \
  LALRunningMedian.c RngMedBias.c
MAINTAINERCLEANFILES = MatrixOps.c
EXTRA_DIST = MatrixOps.m4 LALMatrixOps.m4
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep

MatrixOps.c: MatrixOps.m4 LALMatrixOps.m4
	cd $(srcdir) && @M4@ MatrixOps.m4 > MatrixOps.c
@DVIMKFILE@
