## Process this file with automake to produce Makefile.in
TEXFILES = VectorOps.tex
DVIFILES = VectorOps.dvi
SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = VectorOps.dvi VectorOps.pdf
EXTRA_DIST = $(TEXFILES)
dvi-local: $(DVIFILES)
.tex.dvi:; @LATEX@ $< ; @LATEX@ $<
