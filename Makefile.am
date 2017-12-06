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

BUILT_SOURCES = $(PACKAGE)rc
sysconf_DATA = $(PACKAGE)rc
MOSTLYCLEANFILES = $(PACKAGE)rc
$(PACKAGE)rc: Makefile
	$(AM_V_GEN)echo "# source this file to access $(PACKAGE_NAME)" >$@; \
	for dir in $(SUBDIRS); do \
		echo "eval \`sh -c \"case \$$0 in *csh) echo source $(sysconfdir)/$${dir}-user-env.csh;; *) echo . $(sysconfdir)/$${dir}-user-env.sh;; esac\"\`" >>$@; \
	done

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

.PHONY: print-build-success
all-local: print-build-success
print-build-success:
	@echo "";\
	echo "==================================================";\
	echo "$(PACKAGE_NAME) has now been successfully built.";\
	echo "Now run 'make install' to install $(PACKAGE_NAME).";\
	echo "==================================================";\
	echo ""
.PHONY: print-install-success
install-exec-local: print-install-success
print-install-success:
	@echo "";\
	echo "==================================================";\
	echo "$(PACKAGE_NAME) has now been successfully installed.";\
	echo "To set up your environment, please add the line";\
	echo "  . $(sysconfdir)/$(PACKAGE)rc";\
	echo "to ~/.profile for Bourne shells (e.g. bash), or";\
	echo "  source $(sysconfdir)/$(PACKAGE)rc";\
	echo "to ~/.login for C shells (e.g. tcsh).";\
	echo "==================================================";\
	echo ""
