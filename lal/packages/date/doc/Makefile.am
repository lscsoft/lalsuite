## Process this file with automake to produce Makefile.in
TEXFILES = Date.tex
DVIFILES = Date.dvi
#SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = *.log *.aux
CLEANFILES = Date.dvi Date.pdf
EXTRA_DIST = $(TEXFILES) DateH.tex DateStringC.tex JulianC.tex LMST1C.tex \
             SecsToLALDateC.tex UtoGPSC.tex UtimeC.tex UtoGPSC.tex

dvi-local: $(EXTRA_DIST)
	for file in $(EXTRA_DIST) ; do \
          if [ ! -f $$file ] ; then cp -f $(srcdir)/$$file $$file ; fi ; \
        done
	for texfile in $(TEXFILES) ; do \
          @LATEX@ $$texfile ; @LATEX@ $$texfile ; \
        done

#$(DVIFILES): $(TEXFILES)
#dvi-local: $(DVIFILES)
#.tex.dvi:; @LATEX@ $< ; @LATEX@ $<
