## Process this file with automake to produce Makefile.in

SUBDIRS = laldoc autodoc

man_MANS = lal.3 lal-config.1

INSTALLFILES = lalspec.dvi lalspec.ps lalspec.pdf lsd.dvi lsd.ps lsd.pdf
DISTCLEANFILES = $(INSTALLFILES)
MAINTAINERCLEANFILES = lal-install.tex lal-preinstall.tex \
		       lal-install.txt lal-preinstall.txt \
		       lal-install.sh lal-preinstall.sh
BUILT_SOURCES = lal-install.tex lal-preinstall.tex \
		lal-install.txt lal-preinstall.txt \
		lal-install.sh lal-preinstall.sh
EXTRA_DIST = $(man_MANS) ligodcc.cls lalspec.tex lsd.tex lal.sty \
             lsdFigDirStructure.eps merger.eps \
             lsdFigDirStructure.pdf merger.pdf \
	     lal-install.tex lal-preinstall.tex \
	     lal-install.txt lal-preinstall.txt \
	     lal-install.sh lal-preinstall.sh \
	     lal-install-x.sh lal-preinstall-x.sh \
	     scrstrip.sed texstrip.sed txtstrip.sed

lal-preinstall.sh: lal-preinstall-x.sh scrstrip.sed
	cd $(srcdir) && sed -f scrstrip.sed lal-preinstall-x.sh > lal-preinstall.sh
	cd $(srcdir) && chmod +x lal-preinstall.sh

lal-preinstall.tex: lal-preinstall-x.sh texstrip.sed
	cd $(srcdir) && sed -f texstrip.sed lal-preinstall-x.sh > lal-preinstall.tex

lal-preinstall.txt: lal-preinstall-x.sh txtstrip.sed
	cd $(srcdir) && sed -f txtstrip.sed lal-preinstall-x.sh > lal-preinstall.txt

lal-install.sh: lal-install-x.sh scrstrip.sed
	cd $(srcdir) && sed -f scrstrip.sed lal-install-x.sh > lal-install.sh

lal-install.tex: lal-install-x.sh texstrip.sed
	cd $(srcdir) && sed -f texstrip.sed lal-install-x.sh > lal-install.tex

lal-install.txt: lal-install-x.sh txtstrip.sed
	cd $(srcdir) && sed -f txtstrip.sed lal-install-x.sh > lal-install.txt

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
