## Process this file with automake to produce Makefile.in
TEXFILES = Overlap.tex Dirichlet.tex
DVIFILES = Overlap.dvi Dirichlet.dvi
SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = *.dvi
EXTRA_DIST = $(TEXFILES) *.ps
dvi-local: $(DVIFILES)
.tex.dvi:; @LATEX@ $< ; @LATEX@ $<
