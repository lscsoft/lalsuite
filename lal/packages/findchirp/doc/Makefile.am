## Process this file with automake to produce Makefile.in
TEXFILES = FindChirpExch.tex
DVIFILES = FindChirpExch.dvi
SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = FindChirpExch.dvi FindChirpExch.pdf
EXTRA_DIST = $(TEXFILES)
dvi-local: $(DVIFILES)
.tex.dvi:; @LATEX@ $< ; @LATEX@ $<

