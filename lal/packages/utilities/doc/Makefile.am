## Process this file with automake to produce Makefile.in
PACKAGE = utilities
MOSTLYCLEANFILES = main.log main.ilg main.blg main.out
CLEANFILES = main.aux main.bbl main.glo main.ind main.idx \
             main.lof main.lot main.toc
DISTCLEANFILES = main.dvi main.ps main.pdf
EXTRA_DIST = main.tex $(srcdir)/$(PACKAGE)*

dvi-local:
	$(top_builddir)/doc/laldoc/ldoc-pkg.sh $(PACKAGE) $(srcdir) $(top_builddir) $(top_srcdir)
