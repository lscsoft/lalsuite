## Process this file with automake to produce Makefile.in
TEXFILES = Window.tex
DVIFILES = Window.dvi
SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = Window.dvi Window.pdf
EXTRA_DIST = $(TEXFILES) WindowFig.eps WindowFig.pdf
dvi-local: $(DVIFILES)
.tex.dvi:; @LATEX@ $< ; @LATEX@ $<
