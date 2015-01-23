if DOXYGEN

if HAVE_GIT_REPO
doxygen_include_gitid = true
endif # HAVE_GIT_REPO

doxygen_quiet   = $(doxygen_quiet_$(V))
doxygen_quiet_  = $(doxygen_quiet_$(AM_DEFAULT_VERBOSITY))
doxygen_quiet_0 = YES
doxygen_quiet_1 = NO

html-local: Makefile $(CONFIG_CLEAN_FILES)
	$(AM_V_at)export doxygen_version="$(VERSION)"; \
	if test "x$(doxygen_include_gitid)" = xtrue; then \
		gitid=`cd $(abs_top_srcdir) && $(GIT) log -1 --pretty='format:%h'`; \
		if test "x$${gitid}" != x; then \
			unclean=`cd $(abs_top_srcdir) && ( $(GIT) diff --quiet HEAD || printf %s -UNCLEAN )`; \
			doxygen_version="$${doxygen_version}-$${gitid}$${unclean}"; \
		fi; \
	fi; \
	export doxygen_tagfiles="$(DOXYGEN_TAGFILES)"; \
	export doxygen_quiet="$(doxygen_quiet)"; \
	export doxygen_input=; \
	for dir in packages lib src test python; do \
		test -d "$(abs_top_srcdir)/$${dir}" && doxygen_input="$${doxygen_input} $(abs_top_srcdir)/$${dir}"; \
	done; \
	for file in swig/SWIGCommon.i swig/SWIG$(PACKAGE_NAME)Alpha.i swig/SWIG$(PACKAGE_NAME)Omega.i; do \
		test -r "$(abs_top_srcdir)/$${file}" && doxygen_input="$${doxygen_input} $(abs_top_srcdir)/$${file}"; \
	done; \
	$(MKDIR_P) out/ || exit 1; \
	$(DOXYGEN) doxygen.cfg 2>&1; \
	$(PYTHON) $(abs_srcdir)/check_tags.py $(abs_builddir)/out/$(PACKAGE).tag >> warnings.log || exit 1; \
	$(SED) -e '/^$$/d' -e '/^citelist/d' -e '/warning: Duplicate anchor/d' -e '/warning: explicit link request/d' warnings.log > errors.log; \
	rm -f warnings.log || exit 1; \
	cat errors.log >&2 || exit 1; \
	test `cat errors.log | wc -l` -eq 0 || exit 1; \
	rm -f errors.log || exit 1; \
	echo "=== Doxygen documentation build was successful!"; \
	echo "=== $(PACKAGE_NAME) main page is at file://$(abs_builddir)/out/index.html"

install-data-local: install-html-local
install-html-local: html-local
	$(AM_V_at)for dir in `cd out/ && find . -type d`; do \
		echo " $(MKDIR_P) '$(DESTDIR)$(htmldir)/$${dir}'"; \
		$(MKDIR_P) "$(DESTDIR)$(htmldir)/$${dir}" || exit 1; \
		$(PYTHON) $(abs_srcdir)/install_regex.py "$(DESTDIR)$(htmldir)/$${dir}" "$(DOXYGEN_INSTALL_DIRMAP)" > $@.regex || exit 1; \
		for file in `cd out/ && echo $${dir}/*`; do \
			if test -f "out/$${file}"; then \
				echo " $(INSTALL_DATA) out/$${file} '$(DESTDIR)$(htmldir)/$${file}'"; \
				$(SED) -n -f $@.regex "out/$${file}" > $@.tmp || exit 1; \
				$(INSTALL_DATA) $@.tmp "$(DESTDIR)$(htmldir)/$${file}" || exit 1; \
			fi; \
		done; \
	done; \
	rm -f $@.regex $@.tmp

uninstall-local:
	-rm -rf "$(DESTDIR)$(htmldir)/"

clean-local:
	-rm -rf out/
	-rm -f *.log

EXTRA_DIST = \
	check_tags.py \
	filter_i.sed \
	filter_py.py \
	install_regex.py \
	layout.xml \
	mathjax.js \
	src \
	stylesheet.css \
	toplevel.dox \
	$(END_OF_LIST)

endif # DOXYGEN
