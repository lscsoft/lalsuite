## Process this file with automake to produce Makefile.in
TEXFILES = FrameData.tex SpecBuffer.tex DataBuffer.tex
DVIFILES = FrameData.dvi SpecBuffer.dvi DataBuffer.dvi
SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = FrameData.dvi FrameData.pdf \
             SpecBuffer.dvi SpecBuffer.pdf \
             DataBuffer.dvi DataBuffer.pdf
EXTRA_DIST = $(TEXFILES)
dvi-local: $(DVIFILES)
.tex.dvi:; @LATEX@ $< ; @LATEX@ $<

