## Process this file with automake to produce Makefile.in
TEXFILES = LALConstants.tex LALDatatypes.tex LALStatusMacros.tex LALStdlib.tex
DVIFILES = LALConstants.dvi LALDatatypes.dvi LALStatusMacros.dvi LALStdlib.dvi
SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = *.dvi
EXTRA_DIST = $(TEXFILES)
dvi-local: $(DVIFILES)
.tex.dvi:; @LATEX@ $< ; @LATEX@ $<

