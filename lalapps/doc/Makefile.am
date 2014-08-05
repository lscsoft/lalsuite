## Add additional tex sources (texfiles, eps&pdf figures, etc.) here
TEXSOURCEFILES = lalapps.tex lalapps.bib lalapps.bst \
                 animate.tex animate.ps animate.pdf \
                 BankEfficiency.tex BankNumber.tex \
		 findchirp.tex inspinj.tex inspiral-module.tex \
		 inspiral.tex \
                 stochastic.tex olapredfcn.tex \
                 pipeline-module.tex ring.tex

## Don't edit below here (unless you know what you're doing)
INSTALLFILES = lalapps.pdf lalapps.dvi lalapps.ps
MOSTLYCLEANFILES = *.log *.ilg *.blg *.out
CLEANFILES = *.aux *.bbl *.dvi *.glo *.idx *.ind *.lof *.lot *.toc AUTHORS.tex
DISTCLEANFILES = $(INSTALLFILES)
EXTRA_DIST = $(TEXSOURCEFILES)
LATEX = @LATEX@
MKIND = @MKIND@
DVIPS = @DVIPS@
BIBTEX = @BIBTEX@

pkgdocdir=${prefix}/doc/@PACKAGE@-@VERSION@

clean-local: mostlyclean
	-test -z "$(CLEANFILES)" || rm -f $(CLEANFILES)
	-test "`cd $(srcdir) && pwd`" = "`pwd`" || rm -f $(TEXSOURCEFILES)

dvi-local:
	@for file in $(TEXSOURCEFILES) ; do \
	  test -f $$file || @LN_S@ $(srcdir)/$$file . ; \
	done
	test -f AUTHORS.tex || @LN_S@ $(top_srcdir)/AUTHORS AUTHORS.tex
	${LATEX} lalapps
	${BIBTEX} lalapps
	${LATEX} lalapps
	${MKIND} lalapps
	${LATEX} lalapps
	test ! -f lalapps.dvi || ${DVIPS} lalapps -o

install-data-local:
	$(mkinstalldirs) $(DESTDIR)$(pkgdocdir)
	@for file in $(INSTALLFILES) ; do \
	  if test -f $$file ; then \
	    $(INSTALL_DATA) $$file $(DESTDIR)$(pkgdocdir)/$$file ; \
	    echo "$(INSTALL_DATA) $$file $(DESTDIR)$(pkgdocdir)/$$file" ; \
	  fi; \
	done

uninstall-local:
	@for file in $(INSTALLFILES) ; do \
          if test -f $(DESTDIR)$(pkgdocdir)/$$file ; then \
            rm -f $(DESTDIR)$(pkgdocdir)/$$file; \
            echo "rm -f $(DESTDIR)$(pkgdocdir)/$$file"; \
          fi; \
        done
