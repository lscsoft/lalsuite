## Process this file with automake to produce Makefile.in
TEXFILES = Comm.tex
DVIFILES = Comm.dvi
SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = Comm.dvi Comm.pdf
EXTRA_DIST = $(TEXFILES)
dvi-local: $(DVIFILES)
.tex.dvi:; @LATEX@ $< ; @LATEX@ $<

