## Process this file with automake to produce Makefile.in
TEXFILES = Random.tex FindRoot.tex Sort.tex
DVIFILES = Random.dvi FindRoot.dvi Sort.dvi
SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = *.dvi
EXTRA_DIST = $(TEXFILES) HeapSortC.tex HeapSortD.tex SortH.tex SortTestC.tex
dvi-local: $(DVIFILES)
.tex.dvi:; @LATEX@ $< ; @LATEX@ $<

