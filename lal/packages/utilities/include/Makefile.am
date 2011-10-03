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
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex

all-local: include-link

dvi-local: .dvi-dep

@INCLUDEMKFILE@
@DVIMKFILE@

swig_headers = $(pkginclude_HEADERS)
CLEANFILES += $(swig_clean)
@SWIG_HEADER_MK@
