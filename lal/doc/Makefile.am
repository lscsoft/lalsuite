## Process this file with automake to produce Makefile.in

SUBDIRS = laldoc autodoc

MOSTLYCLEANFILES = lsd.log lsd.ilg lsd.blg lsd.out\
                   lalspec.log lalspec.ilg lalspec.blg lalspec.out
CLEANFILES = lsd.aux lsd.bbl lsd.glo lsd.ind lsd.idx lsd.lof lsd.lot lsd.toc\
             lalspec.aux lalspec.bbl lalspec.glo lalspec.ind lalspec.idx\
             lalspec.lof lalspec.lot lalspec.toc
DISTCLEANFILES = lalspec.dvi lalspec.ps lalspec.pdf lsd.dvi lsd.ps lsd.pdf
EXTRA_DIST = ligodcc.cls lalspec.tex lsd.tex lal.sty \
             lsdFigDirStructure.eps schedule.eps \
             lsdFigDirStructure.pdf schedule.pdf

pkgdocdir=${prefix}/doc/@PACKAGE@-@VERSION@

dvi-local:
	@LATEX@ lalspec && @LATEX@ lalspec && @LATEX@ lalspec

install-data-local:
	$(mkinstalldirs) $(DESTDIR)$(pkgdocdir)
	if test -f lsd.pdf; then \
          $(INSTALL_DATA) lsd.pdf $(DESTDIR)$(pkgdocdir)/lsd.pdf ; \
        fi
	if test -f lsd.dvi; then \
	  $(INSTALL_DATA) lsd.dvi $(DESTDIR)$(pkgdocdir)/lsd.dvi ; \
        fi
	if test -f lalspec.pdf; then \
          $(INSTALL_DATA) lalspec.pdf $(DESTDIR)$(pkgdocdir)/lalspec.pdf ; \
        fi
	if test -f lalspec.dvi; then \
          $(INSTALL_DATA) lalspec.dvi $(DESTDIR)$(pkgdocdir)/lalspec.dvi ; \
        fi

uninstall-local:
	rm -f $(DESTDIR)$(pkgdocdir)/lsd.pdf
	rm -f $(DESTDIR)$(pkgdocdir)/lsd.dvi
	rm -f $(DESTDIR)$(pkgdocdir)/lalspec.pdf
	rm -f $(DESTDIR)$(pkgdocdir)/lalspec.dvi
