## Process this file with automake to produce Makefile.in
TEXFILES = LALError.tex LALMalloc.tex PrintVector.tex
DVIFILES = LALError.dvi LALMalloc.dvi PrintVector.dvi
SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = *.dvi
EXTRA_DIST = $(TEXFILES)
dvi-local: $(DVIFILES)
.tex.dvi:; @LATEX@ $< ; @LATEX@ $<

