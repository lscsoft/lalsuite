## Process this file with automake to produce Makefile.in
TEXFILES = ComplexFFT.tex RealFFT.tex
DVIFILES = ComplexFFT.dvi RealFFT.dvi
SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = ComplexFFT.dvi ComplexFFT.pdf RealFFT.dvi RealFFT.pdf
EXTRA_DIST = $(TEXFILES)
dvi-local: $(DVIFILES)
.tex.dvi:; @LATEX@ $< ; @LATEX@ $<
