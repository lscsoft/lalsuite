pkginclude_HEADERS = \
	BitField.h \
	CoarseGrainFrequencySeries.h \
	Dirichlet.h \
	FindRoot.h \
	GSLSupport.h \
	Integrate.h \
	Interpolate.h \
	LALAdaptiveRungeKutta4.h \
	LALRunningMedian.h \
	MatrixUtils.h \
	ODE.h \
	Random.h \
	RngMedBias.h \
	Sort.h \
	SphericalHarmonics.h

CLEANFILES =
MOSTLYCLEANFILES =
DISTCLEANFILES =

all-local: include-link

@INCLUDEMKFILE@

swig_headers = $(pkginclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
