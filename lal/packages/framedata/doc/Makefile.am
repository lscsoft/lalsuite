## Process this file with automake to produce Makefile.in
TEXFILES = FrameData.tex SpecBuffer.tex DataBuffer.tex
DVIFILES = FrameData.dvi SpecBuffer.dvi DataBuffer.dvi
SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = *.dvi
EXTRA_DIST = $(TEXFILES)
dvi-local: $(DVIFILES)
.tex.dvi:; @LATEX@ $< ; @LATEX@ $<

