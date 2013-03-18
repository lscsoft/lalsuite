ACLOCAL_AMFLAGS = -I gnuscripts

# Override automake's default rule for rebuilding aclocal.m4, to ensure that
# ./configure is successfully recreated after certain build system changes which
# require other autotools (e.g. autoheader) to be re-run. See Redmine issue #728.
$(ACLOCAL_M4): $(am__aclocal_m4_deps)
	$(am__cd) $(srcdir) && $(SHELL) ./00boot

# Override automake's default rule for rebuilding ./config.status, to ensure that
# all changes to the top-level configuration (especially to environment variables it
# exports) are propagated to the library-level configurations. See Redmine issue #728.
$(top_builddir)/config.status: $(top_srcdir)/configure $(CONFIG_STATUS_DEPENDENCIES)
	$(SHELL) $(top_srcdir)/configure $(ac_configure_args)

SUBDIRS = \
	doxygen \
	gnuscripts \
	@subdirs@

distcleancheck_listfiles = find . -type f \! -name "LAL*VCSInfo.h" -print

vcs-clean: maintainer-clean
	rm -f Makefile.in
	rm -f configure
	rm -f aclocal.m4
	rm -f config.cache
	rm -f config.status
	rm -f config.log
	rm -f gnuscripts/Makefile.in

cvs-clean: vcs-clean
