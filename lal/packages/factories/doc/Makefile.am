## Process this file with automake to produce Makefile.in
PACKAGE = factories
MOSTLYCLEANFILES = main.log main.ilg main.blg main.out
CLEANFILES = main.aux main.bbl main.glo main.ind main.idx \
             main.lof main.lot main.toc
DISTCLEANFILES = main.dvi main.ps main.pdf
EXTRA_DIST = main.tex $(PACKAGE)*
ADOCDIR = $(top_srcdir)/doc/autodoc
LDOCDIR = $(top_srcdir)/doc/laldoc
LALDOC = ls *.[ch] | xargs -i $(LDOCDIR)/laldoc {} logfile $(ADOCDIR)/ ; mv *.tex $(ADOCDIR)

dvi-local:
	cd $(srcdir)/../include; $(LALDOC)
	cd $(srcdir)/../src; $(LALDOC)
	cd $(srcdir)/../test; $(LALDOC)
	cd $(ADOCDIR) ; rm -f $(PACKAGE)*
	cd $(srcdir) ; @LN_S@ `pwd`/$(PACKAGE)* $(ADOCDIR)
	cp main.* $(top_srcdir)/doc/autodoc
	cd $(ADOCDIR); @LATEX@ main.tex; @MKIND@ main.idx; @LATEX@ main.tex; @LATEX@ main.tex
	rm -f $(ADOCDIR)/main.tex
	mv $(ADOCDIR)/main.* $(srcdir)
