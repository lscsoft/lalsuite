## Process this file with automake to produce Makefile.in
TEXFILES = Thresholds.tex TFTransform.tex ExcessPower.tex
DVIFILES = Thresholds.dvi TFTransform.dvi ExcessPower.dvi
SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = Thresholds.dvi TFTransform.dvi ExcessPower.dvi \
             Thresholds.pdf TFTransform.pdf ExcessPower.pdf
EXTRA_DIST = $(TEXFILES)
dvi-local: $(DVIFILES)
.tex.dvi:; @LATEX@ $< ; @LATEX@ $<
