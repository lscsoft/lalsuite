## Process this file with automake to produce Makefile.in

SUBDIRS = laldoc autodoc

MOSTLYCLEANFILES = lsd.log lsd.ilg lsd.blg lsd.out
CLEANFILES = lsd.aux lsd.bbl lsd.glo lsd.ind lsd.idx lsd.lof lsd.lot lsd.toc
DISTCLEANFILES = lsd.dvi lsd.ps lsd.pdf
EXTRA_DIST = lsd.tex lal.sty \
             lsdFigDirStructure.eps schedule.eps \
             lsdFigDirStructure.pdf schedule.pdf

pkgdocdir=@prefix@/doc/@PACKAGE@-@VERSION@

install-data-local:
	$(mkinstalldirs) $(DESTDIR)$(pkgdocdir)
	if test -f lsd.pdf; then \
          $(INSTALL_DATA) lsd.pdf $(DESTDIR)$(pkgdocdir)/lsd.pdf ; \
        fi
	if test -f lsd.dvi; then \
	  $(INSTALL_DATA) lsd.dvi $(DESTDIR)$(pkgdocdir)/lsd.dvi ; \
        fi

uninstall-local:
	rm -f $(DESTDIR)$(pkgdocdir)/lsd.pdf
	rm -f $(DESTDIR)$(pkgdocdir)/lsd.dvi
