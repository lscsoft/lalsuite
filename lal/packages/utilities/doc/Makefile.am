## Process this file with automake to produce Makefile.in
TEXFILES = Random.tex FindRoot.tex Integrate.tex Interpolate.tex Sort.tex
DVIFILES = Random.dvi FindRoot.dvi Integrate.dvi Interpolate.dvi Sort.dvi
#SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = Random.dvi FindRoot.dvi Integrate.dvi Interpolate.dvi Sort.dvi \
             Random.pdf FindRoot.pdf Integrate.pdf Interpolate.pdf Sort.pdf
EXTRA_DIST = $(TEXFILES) HeapSortC.tex HeapSortD.tex SortH.tex SortTestC.tex


dvi-local: $(EXTRA_DIST)
	for file in $(EXTRA_DIST) ; do \
          if [ ! -f $$file ] ; then cp -f $(srcdir)/$$file $$file ; fi ; \
        done
	for texfile in $(TEXFILES) ; do \
          @LATEX@ $$texfile ; @LATEX@ $$texfile ; \
        done

#dvi-local: $(DVIFILES)
#.tex.dvi:; @LATEX@ $< ; @LATEX@ $<
