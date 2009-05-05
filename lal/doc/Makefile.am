## Process this file with automake to produce Makefile.in

SUBDIRS = laldoc autodoc

man_MANS = lal.3 lal-config.1

INSTALLFILES = lalspec.dvi lalspec.ps lalspec.pdf lsd.dvi lsd.ps lsd.pdf
DISTCLEANFILES = $(INSTALLFILES)
EXTRA_DIST = $(man_MANS) ligodcc.cls lalspec.tex lsd.tex lal.sty \
             lsdFigDirStructure.eps merger.eps \
						 lsdFigDirStructure.pdf merger.pdf

install-data-local:
	$(mkinstalldirs) $(DESTDIR)$(pkgdocdir)
	@for file in $(INSTALLFILES) ; do \
          if test -f $$file ; then \
            $(INSTALL_DATA) $$file $(DESTDIR)$(pkgdocdir)/$$file; \
            echo "$(INSTALL_DATA) $$file $(DESTDIR)$(pkgdocdir)/$$file"; \
          fi; \
        done

uninstall-local:
	@for file in $(INSTALLFILES) ; do \
          if test -f $(DESTDIR)$(pkgdocdir)/$$file ; then \
            rm -f $(DESTDIR)$(pkgdocdir)/$$file; \
	    echo "rm -f $(DESTDIR)$(pkgdocdir)/$$file"; \
          fi; \
        done
