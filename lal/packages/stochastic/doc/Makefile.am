## Process this file with automake to produce Makefile.in
DOCSOURCEFILES = stochastic.tex \
  stochasticDirichletFig1.eps stochasticDirichletFig1.pdf \
  stochasticDirichletFig2.eps stochasticDirichletFig2.pdf \
  stochasticDirichletFig3.eps stochasticDirichletFig3.pdf \
  stochasticOmegaGWQuadratic.eps stochasticOmegaGWQuadratic.pdf \
  stochasticOverlapLHOLLO.eps stochasticOverlapLHOLLO.pdf
MOSTLYCLEANFILES = .dvi-dep
CLEANFILES = main.dvi main.pdf main.ps
EXTRA_DIST = $(DOCSOURCEFILES)
dvi-local: .dvi-dep
distclean-local:
	-rm -fr .adoc
