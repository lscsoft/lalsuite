## Process this file with automake to produce Makefile.in
TEXFILES = Date.tex DateH.tex DateStringC.tex JulianC.tex LMST1C.tex \
           SecsToLALDateC.tex UtoGPSC.tex UtimeC.tex UtoGPSC.tex
DVIFILES = Date.dvi
SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = Date.dvi Date.pdf
EXTRA_DIST = $(TEXFILES)
$(DVIFILES): $(TEXFILES)
dvi-local: $(DVIFILES)
.tex.dvi:; @LATEX@ $< ; @LATEX@ $<
