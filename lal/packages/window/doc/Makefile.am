## Process this file with automake to produce Makefile.in
TEXFILES = Window.tex
DVIFILES = Window.dvi
SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = Window.dvi Window.pdf
EXTRA_DIST = $(TEXFILES) WindowFig.eps WindowFig.pdf
dvi-local: $(EXTRA_DIST)
	for file in $(EXTRA_DIST) ; do \
          if [ ! -f $$file ] ; then cp -f $(srcdir)/$$file $$file ; fi ; \
        done
	for texfile in $(TEXFILES) ; do \
          @LATEX@ $$texfile ; @LATEX@ $$texfile ; \
        done
#dvi-local: $(DVIFILES)
#.tex.dvi:; @LATEX@ $< ; @LATEX@ $<
