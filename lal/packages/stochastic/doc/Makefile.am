## Process this file with automake to produce Makefile.in
DOCSOURCEFILES = stochastic.tex \
  stochasticOmegaGWQuadratic.eps stochasticOmegaGWQuadratic.pdf \
  stochasticOverlapGEO600TAMA300.eps stochasticOverlapGEO600TAMA300.pdf \
  stochasticOverlapLHOGEO600.eps stochasticOverlapLHOGEO600.pdf \
  stochasticOverlapLHOLLO.eps stochasticOverlapLHOLLO.pdf \
  stochasticOverlapLHOTAMA300.eps stochasticOverlapLHOTAMA300.pdf \
  stochasticOverlapLHOVIRGO.eps stochasticOverlapLHOVIRGO.pdf \
  stochasticOverlapLLOGEO600.eps stochasticOverlapLLOGEO600.pdf \
  stochasticOverlapLLOTAMA300.eps stochasticOverlapLLOTAMA300.pdf \
  stochasticOverlapLLOVIRGO.eps stochasticOverlapLLOVIRGO.pdf \
  stochasticOverlapVIRGOGEO600.eps stochasticOverlapVIRGOGEO600.pdf \
  stochasticOverlapVIRGOTAMA300.eps stochasticOverlapVIRGOTAMA300.pdf \
  omega_binary.eps omega_binary.pdf omega_pulsar.eps omega_pulsar.pdf
MOSTLYCLEANFILES = .dvi-dep
CLEANFILES = main.dvi main.pdf main.ps
EXTRA_DIST = $(DOCSOURCEFILES)
dvi-local: .dvi-dep
distclean-local:
	-rm -fr .adoc
@DOCMKFILE@
