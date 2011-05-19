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

noinst_HEADERS = MatrixOps_source.c

MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex

dvi-local: .dvi-dep

@DVIMKFILE@
