## Process this file with automake to produce Makefile.in
PACKAGE = date

TEXFILES = Date.tex
DVIFILES = Date.dvi
#SUFFIXES = .tex .dvi
MOSTLYCLEANFILES = main.log main.ilg main.blg main.out \
                   *.log *.aux
CLEANFILES = main.aux main.bbl main.glo main.ind main.idx \
             main.lof main.lot main.toc \
             Date.dvi Date.pdf
DISTCLEANFILES = main.dvi main.ps main.pdf
MORE_EXTRA_DIST = $(TEXFILES) DateH.tex DateStringC.tex JulianC.tex LMST1C.tex \
                  SecsToLALDateC.tex UtoGPSC.tex UtimeC.tex UtoGPSC.tex
EXTRA_DIST = main.tex $(srcdir)/$(PACKAGE)* \
             $(MORE_EXTRA_DIST)

dvi-local:
	$(top_builddir)/doc/laldoc/ldoc-pkg.sh $(PACKAGE) $(srcdir) $(top_builddir) $(top_srcdir)
	for file in $(MORE_EXTRA_DIST) ; do \
          if [ ! -f $$file ] ; then cp -f $(srcdir)/$$file $$file ; fi ; \
        done
	for texfile in $(TEXFILES) ; do \
          @LATEX@ $$texfile ; @LATEX@ $$texfile ; \
        done

#$(DVIFILES): $(TEXFILES)
#dvi-local: $(DVIFILES)
#.tex.dvi:; @LATEX@ $< ; @LATEX@ $<
