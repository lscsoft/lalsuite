## Process this file with automake to produce Makefile.in
EXTRA_DIST = doc.mk dvi.mk include.mk main.tex.in \
             lal.m4 lallibtool.m4 lalam15patch.m4 pkg.m4 old-pkg.m4 \
             acconfig.h lal-user-env.sh lal-user-env.csh \
	     capture-uninstall.sh

CLEANFILES = uninstall-lal.sh

sysconf_DATA = lal-user-env.sh lal-user-env.csh

install-exec-local: uninstall-lal-sh-target
	$(mkinstalldirs) $(DESTDIR)$(sbindir)
	$(INSTALL_SCRIPT) uninstall-lal.sh $(DESTDIR)$(sbindir)/uninstall-lal.sh

uninstall-lal-sh-target:
	top_builddir=$(top_builddir) ; export top_builddir ; \
        $(SHELL) $(srcdir)/capture-uninstall.sh uninstall-lal.sh
