## Process this file with automake to produce Makefile.in
DOCSOURCEFILES = utilities.tex \
  utilitiesDirichletFig1.eps utilitiesDirichletFig1.pdf \
  utilitiesDirichletFig2.eps utilitiesDirichletFig2.pdf \
  utilitiesDirichletFig3.eps utilitiesDirichletFig3.pdf
MOSTLYCLEANFILES = .dvi-dep
CLEANFILES = main.dvi main.pdf main.ps
EXTRA_DIST = $(DOCSOURCEFILES)
dvi-local: .dvi-dep
distclean-local:
	-rm -fr .adoc
@DOCMKFILE@
