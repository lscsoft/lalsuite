## Process this file with automake to produce Makefile.in
TEXFILES = Random.tex FindRoot.tex
DVIFILES = Random.dvi FindRoot.dvi
SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = *.dvi
EXTRA_DIST = $(TEXFILES)
dvi-local: $(DVIFILES)
.tex.dvi:; @LATEX@ $< ; @LATEX@ $<

