## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = Random.h FindRoot.h Integrate.h Interpolate.h \
  Sort.h ODE.h Dirichlet.h CoarseGrainFrequencySeries.h MatrixUtils.h \
  LALRunningMedian.h RngMedBias.h SFTClean.h NormalizeSFTRngMed.h
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
all-local: include-link
dvi-local: .dvi-dep
@INCLUDEMKFILE@
@DVIMKFILE@
