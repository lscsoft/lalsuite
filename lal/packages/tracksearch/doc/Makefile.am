## Process this file with automake to produce Makefile.in
TEXFILES = TrackSearch.tex
DVIFILES = TrackSearch.dvi
SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = TrackSearch.dvi TrackSearch.pdf
EXTRA_DIST = $(TEXFILES)
dvi-local: $(DVIFILES)
.tex.dvi:; @LATEX@ $< ; @LATEX@ $<

@DOCMKFILE@
