INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

noinst_LTLIBRARIES = libutilities.la
libutilities_la_SOURCES = \
	CoarseGrainFrequencySeries.c \
	DetInverse.c \
	DetInverseInternal.c \
	Dirichlet.c \
	Eigen.c \
	EigenInternal.c \
	FindRoot.c \
	GSLSupport.c \
	HeapSort.c \
	Integrate.c \
	Interpolate.c \
	LALRunningMedian.c \
	MatrixOps.c \
	MersenneRandom.c \
	ODE.c \
	Random.c \
	RngMedBias.c

MAINTAINERCLEANFILES = MatrixOps.c
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex

EXTRA_DIST = \
	LALMatrixOps.m4 \
	MatrixOps.m4

dvi-local: .dvi-dep

MatrixOps.c: MatrixOps.m4 LALMatrixOps.m4
	cd $(srcdir) && @M4@ MatrixOps.m4 > MatrixOps.c

@DVIMKFILE@
