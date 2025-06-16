BUILT_SOURCES =
sysconf_DATA =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_top.am

EXTRA_DIST += CONTRIBUTING.md

ACLOCAL_AMFLAGS = -I gnuscripts

SUBDIRS = @subdirs@ doxygen wheel

.PHONY: wheel
wheel:
	$(MAKE) -C wheel wheel

# Override automake's default rule for rebuilding aclocal.m4, to ensure that
# ./configure is successfully recreated after certain build system changes which
# require other autotools (e.g. autoheader) to be re-run. See Redmine issue #728.
$(lalsuite_package_configure_deps):
$(ACLOCAL_M4): $(am__aclocal_m4_deps) $(lalsuite_package_configure_deps)
	$(am__cd) $(srcdir) && $(SHELL) ./00boot

# Override automake's default rule for rebuilding ./config.status, to ensure that
# all changes to the top-level configuration (especially to environment variables it
# exports) are propagated to the library-level configurations. See Redmine issue #728.
$(top_builddir)/config.status: $(top_srcdir)/configure $(CONFIG_STATUS_DEPENDENCIES)
	$(SHELL) $(top_srcdir)/configure $(ac_configure_args)

BUILT_SOURCES += AUTHORS
MOSTLYCLEANFILES += AUTHORS
AUTHORS: $(lalsuite_package_AUTHORS_deps)
	$(AM_V_GEN)cat $(lalsuite_package_AUTHORS_deps) | LC_ALL=C sort | uniq > $@.tmp && mv -f $@.tmp $@

.PHONY: cscopelist-subdirs
cscopelist-am: cscopelist-subdirs
cscopelist-subdirs:
	for subdir in @subdirs@; do \
		sed "s|^|$${subdir}/|" $${subdir}/cscope.files >> $(top_builddir)/cscope.files; \
	done

user_environment = \
	`for dir in @subdirs@; do printf 'source %s; ' "$(sysconfdir)/$${dir}-user-env"; done` \
	$(END_OF_LIST)

if !SYSCONFDIR_IS_ETC

BUILT_SOURCES += $(PACKAGE)rc
sysconf_DATA += $(PACKAGE)rc
MOSTLYCLEANFILES += $(PACKAGE)rc
$(PACKAGE)rc: Makefile
	$(AM_V_GEN)echo "# source this file to access $(PACKAGE_NAME) from Bourne or C shells" >$@; \
	echo "expr \"X\$$0\" : '^X.*csh' >/dev/null && source $(sysconfdir)/$(PACKAGE)-user-env.csh || . $(sysconfdir)/$(PACKAGE)-user-env.sh" >>$@

endif

if DOXYGEN

html-local:
	@echo "$(subdir)/Makefile: Doxygen documentation was built successfully!"; \
	echo "$(subdir)/Makefile: $(PACKAGE_NAME) main page is at file://$(abs_builddir)/doxygen/out/index.html"

install-html-local:
	@echo "$(subdir)/Makefile: Doxygen documentation was installed successfully!"; \
	echo "$(subdir)/Makefile: installed $(PACKAGE_NAME) main page is at file://$(DESTDIR)$(htmldir)/index.html"

uninstall-local:
	-rm -rf "$(DESTDIR)$(htmldir)/"

endif
