## Process this file with automake to produce Makefile.in
SUBDIRS = include packages
LIBLAL = $(srcdir)/lib/liblal.a
CLEANFILES = $(LIBLAL) $(srcdir)/include/*.h
EXTRA_DIST = $(srcdir)/lib/README lalspec.pdf

#disable autoheader
AUTOHEADER=echo

install-exec-local:
	@$(NORMAL_INSTALL)
	$(mkinstalldirs) $(DESTDIR)$(pkglibdir)
	echo "$(INSTALL) $(LIBLAL) $(DESTDIR)$(pkglibdir)/liblal.a"
	$(INSTALL) $(LIBLAL) $(DESTDIR)$(pkglibdir)/liblal.a

uninstall-local:
	@$(NORMAL_UNINSTALL)
	rm -f $(DESTDIR)$(pkglibdir)/liblal.a

