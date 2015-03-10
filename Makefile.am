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

BUILT_SOURCES = $(PACKAGE)rc
sysconf_DATA = $(PACKAGE)rc
MOSTLYCLEANFILES = $(PACKAGE)rc
$(PACKAGE)rc: Makefile
	$(AM_V_GEN)echo "# source this file to access $(PACKAGE_NAME)" >$@; \
	for dir in $(SUBDIRS); do \
		echo "source $(sysconfdir)/$${dir}-user-env.\`sh -c \"case \$$0 in *csh) echo csh;; *) echo sh;; esac\"\`" >>$@; \
	done

if DOXYGEN

html-local: Makefile $(CONFIG_CLEAN_FILES)
	$(AM_V_at)echo "=== Doxygen documentation was built successfully!"; \
	echo "=== $(PACKAGE_NAME) main page is at file://$(abs_builddir)/index.html"

install-html-local: html-local
	$(AM_V_at)echo " $(MKDIR_P) '$(DESTDIR)$(htmldir)'"; \
	$(MKDIR_P) "$(DESTDIR)$(htmldir)" || exit 1; \
	echo " $(INSTALL_DATA) index.html '$(DESTDIR)$(htmldir)/index.html'"; \
	$(SED) -e 's|\(lal[a-z]*\)/doxygen/out|../\1|g' index.html > $@.tmp || exit 1; \
	$(INSTALL_DATA) $@.tmp "$(DESTDIR)$(htmldir)/index.html" || exit 1; \
	rm -f $@.tmp; \
	echo "=== Doxygen documentation was installed successfully!"; \
	echo "=== Installed $(PACKAGE_NAME) main page is at file://$(DESTDIR)$(htmldir)/index.html"

uninstall-local:
	-rm -rf "$(DESTDIR)$(htmldir)/"

endif
