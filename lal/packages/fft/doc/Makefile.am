## Process this file with automake to produce Makefile.in
TEXFILES = ComplexFFT.tex RealFFT.tex
DVIFILES = ComplexFFT.dvi RealFFT.dvi
SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = *.dvi
EXTRA_DIST = $(TEXFILES)
dvi-local: $(DVIFILES)
.tex.dvi:; @LATEX@ $< ; @LATEX@ $<
