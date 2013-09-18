DOXY_HTML_FILES = \
	doxy-html-config.cfg \
	doxy-html-mathjax.js \
	doxy-html-layout.xml

DOXY_PDF_FILES = \
	doxy-pdf-config.cfg  \
	doxy-pdf-extra.sty \
	doxy-pdf-layout.xml

EXTRA_DIST = \
	$(DOXY_HTML_FILES) \
	$(DOXY_PDF_FILES)

SUBDIRS = src

.PHONY: html pdf

DOXYINPUT = \
	lal \
	lalburst \
	laldetchar \
	lalframe \
	lalinference \
	lalinspiral \
	lalmetaio \
	lalpulsar \
	lalsimulation \
	lalstochastic \
	lalxml \
	lalapps

override doxy_inputdirs = $(abs_srcdir)/src \
	$(addprefix $(abs_top_srcdir)/,$(DOXYINPUT))

override lalsuite_version = $(shell \
	version="$(PACKAGE_VERSION)"; \
	if test "x$(GIT)" != x; then \
		hash=`cd $(abs_top_srcdir) && $(GIT) rev-parse --quiet HEAD 2>/dev/null`; \
		if test "x$${hash}" != x; then \
			taghash=`cd $(abs_top_srcdir) && $(GIT) rev-parse lalsuite-v${PACKAGE_VERSION}^{commit} 2>/dev/null`; \
			unclean=`cd $(abs_top_srcdir) && ( $(GIT) diff --quiet HEAD || printf %s -UNCLEAN )`; \
			if test "x$${hash}" != "x$${taghash}" || test "x$${unclean}" != x; then \
				version="$${version}-$${hash}$${unclean}"; \
			fi; \
		fi; \
	fi; \
	echo "$${version}" )

override TEXINPUTS = .:$(abs_srcdir):
override pool_size = 4000000

export abs_top_srcdir abs_srcdir abs_builddir
export doxy_inputdirs lalsuite_version TEXINPUTS pool_size

html: $(DOXY_HTML_FILES)
	@echo "===== building doxygen HTML in $(abs_builddir) ======"
	@echo "DOXYINPUT=$(DOXYINPUT)"
	rm -f doxy-html*.log
	if doxygen $(abs_srcdir)/doxy-html-config.cfg; then \
		if test -f html/_formulas.log; then \
			$(LN_S) html/_formulas.log doxy-html-eqns.log; \
		fi; \
	fi
	@echo "===== OK: done building doxygen HTML ====="

pdf: $(DOXY_PDF_FILES)
	@echo "===== building doxygen PDF in $(abs_builddir) ====="
	@echo "DOXYINPUT=$(DOXYINPUT)"
	rm -f refman.pdf doxy-pdf*.log
	if doxygen $(abs_srcdir)/doxy-pdf-config.cfg && ( cd latex && $(MAKE) ); then \
		$(LN_S) latex/refman.pdf refman.pdf; \
		exit 0; \
	else \
		$(LN_S) latex/refman.log doxy-pdf-make.log; \
		exit 1; \
	fi
	@echo "===== OK: done building doxygen PDF ====="

clean-local:
	-rm -fr html/ latex/
	-rm -f  refman.pdf doxy-html*.log doxy-pdf*.log
