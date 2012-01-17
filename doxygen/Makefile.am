DOXY_HTML_FILES = \
	doxy-html-config.cfg \
	doxy-html-aliases.cfg \
	doxy-html-eqns.sty \
	doxy-html-layout.xml

DOXY_PDF_FILES = \
	doxy-pdf-config.cfg  \
	doxy-pdf-aliases.cfg \
	doxy-pdf-header.tex \
	doxy-pdf-layout.xml

EXTRA_DIST = \
	$(DOXY_HTML_FILES) \
	$(DOXY_PDF_FILES)

SUBDIRS = src

.PHONY: html pdf

DOXYINPUT = \
	lal/packages \
	lalburst \
	lalframe \
	lalinference \
	lalinspiral \
	lalmetaio \
	lalpulsar \
	lalsimulation \
	lalstochastic \
	lalxml \
	lalapps/src

override doxy_abs_top_srcdir = $(shell cd $(top_srcdir) && pwd)
override doxy_abs_srcdir = $(shell cd $(srcdir) && pwd)
override doxy_abs_builddir = $(shell cd $(builddir) && pwd)

override doxy_inputdirs = $(doxy_abs_srcdir)/src \
	$(addprefix $(doxy_abs_top_srcdir)/,$(DOXYINPUT))

override lalsuite_version = $(shell \
	cd $(doxy_abs_top_srcdir) || exit 1; \
	lsv=$(PACKAGE_VERSION); \
	if test "x$(GIT)" != x && test -d .git; then \
		lsv="$$lsv-"`$(GIT) log -1 --pretty="format:%h"`; \
		if $(GIT) diff --quiet HEAD; then \
			: ; \
		else \
			lsv="$$lsv-unclean"; \
		fi; \
	fi; \
	echo $$lsv )

override TEXINPUTS = .:$(doxy_abs_srcdir):
override pool_size = 4000000

export doxy_abs_top_srcdir doxy_abs_srcdir doxy_abs_builddir
export doxy_inputdirs lalsuite_version TEXINPUTS pool_size

html: $(DOXY_HTML_FILES)
	@echo "===== building doxygen HTML in $(doxy_abs_builddir) ======"
	@echo "DOXYINPUT=$(DOXYINPUT)"
	rm -f doxy-html*.log
	if doxygen $(doxy_abs_srcdir)/doxy-html-config.cfg; then \
		if test -f html/_formulas.log; then \
			$(LN_S) html/_formulas.log doxy-html-eqns.log; \
		fi; \
		 exit 0; \
	fi
	@echo "===== OK: done building doxygen HTML ====="

pdf: $(DOXY_PDF_FILES)
	@echo "===== building doxygen PDF in $(doxy_abs_builddir) ====="
	@echo "DOXYINPUT=$(DOXYINPUT)"
	rm -f refman.pdf doxy-pdf*.log
	if doxygen $(doxy_abs_srcdir)/doxy-pdf-config.cfg && ( cd latex && $(MAKE) ); then \
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
