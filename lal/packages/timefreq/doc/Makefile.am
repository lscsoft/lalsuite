## Process this file with automake to produce Makefile.in
TEXFILES = Tfr.tex
DVIFILES = Tfr.dvi
SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = Tfr.dvi Tfr.pdf
EXTRA_DIST = $(TEXFILES)
dvi-local: $(DVIFILES)
.tex.dvi:; @LATEX@ $< ; @LATEX@ $<

@DOCMKFILE@
