BUILT_SOURCES =
sysconf_DATA =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_top.am

ACLOCAL_AMFLAGS = -I gnuscripts

SUBDIRS = @subdirs@

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

.PHONY: cscopelist-subdirs
cscopelist-am: cscopelist-subdirs
cscopelist-subdirs:
	for subdir in $(SUBDIRS); do \
		sed "s|^|$${subdir}/|" $${subdir}/cscope.files >> $(top_builddir)/cscope.files; \
	done

user_environment = \
	`for dir in $(SUBDIRS); do printf 'source %s; ' "$(sysconfdir)/$${dir}-user-env"; done` \
	$(END_OF_LIST)

BUILT_SOURCES += $(PACKAGE)rc
sysconf_DATA += $(PACKAGE)rc
MOSTLYCLEANFILES += $(PACKAGE)rc
$(PACKAGE)rc: Makefile
	$(AM_V_GEN)echo "echo '$@ is deprecated; source $(sysconfdir)/$(PACKAGE)-user-env.[c]sh instead'; expr \"X\$$0\" : '^X.*csh' >/dev/null && source $(sysconfdir)/$(PACKAGE)-user-env.csh >/dev/null || . $(sysconfdir)/$(PACKAGE)-user-env.sh >/dev/null" >$@

if DOXYGEN

html-local: Makefile $(CONFIG_CLEAN_FILES)
	$(AM_V_at)set -e; \
	echo "$(subdir)/Makefile: Doxygen documentation was built successfully!"; \
	echo "$(subdir)/Makefile: $(PACKAGE_NAME) main page is at file://$(abs_builddir)/index.html"

install-html-local: html-local
	$(AM_V_at)set -e; \
	echo " $(MKDIR_P) '$(DESTDIR)$(htmldir)'"; \
	$(MKDIR_P) "$(DESTDIR)$(htmldir)"; \
	echo " $(INSTALL_DATA) index.html '$(DESTDIR)$(htmldir)/index.html'"; \
	$(SED) -e 's|\(lal[a-z]*\)/doxygen/out|../\1|g' index.html > $@.tmp; \
	$(INSTALL_DATA) $@.tmp "$(DESTDIR)$(htmldir)/index.html"; \
	rm -f $@.tmp; \
	echo "$(subdir)/Makefile: Doxygen documentation was installed successfully!"; \
	echo "$(subdir)/Makefile: installed $(PACKAGE_NAME) main page is at file://$(DESTDIR)$(htmldir)/index.html"

uninstall-local:
	-rm -rf "$(DESTDIR)$(htmldir)/"

endif
