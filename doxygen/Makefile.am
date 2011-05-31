DOXY_HTML_FILES = doxy-html-config.cfg.in \
		  doxy-html-aliases.cfg \
		  doxy-html-eqns.sty \
		  doxy-html-layout.xml

DOXY_PDF_FILES = doxy-pdf-config.cfg.in  \
		 doxy-pdf-aliases.cfg \
		 doxy-pdf-header.tex \
		 doxy-pdf-layout.xml

EXTRA_DIST = $(DOXY_HTML_FILES) $(DOXY_PDF_FILES)

SUBDIRS = src

.PHONY: html pdf

DOXYINPUT = lal/packages \
            lalframe \
            lalmetaio \
            lalxml \
            lalburst \
            lalinspiral \
            lalpulsar \
            lalstochastic \
            lalapps/src

override doxy_inputdirs = $(abs_srcdir)/src \
                          $(addprefix $(abs_lalsuite_srcdir)/,$(DOXYINPUT))

override lalsuite_version = $(shell \
   cd $(abs_lalsuite_srcdir) || exit 1; \
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

override TEXINPUTS = .:$(abs_srcdir):
override pool_size = 4000000

export doxy_inputdirs lalsuite_version TEXINPUTS pool_size

html: $(DOXY_HTML_FILES)
	@echo "===== building doxygen HTML in $(abs_builddir) ======"
	@echo "DOXYINPUT=$(DOXYINPUT)"
	rm -f doxy-html*.log
	cd $(abs_builddir) && ( \
	   if doxygen doxy-html-config.cfg; then \
	      if test -f html/_formulas.log; then \
	         $(LN_S) html/_formulas.log doxy-html-eqns.log; \
	      fi; \
	      exit 0; \
	   fi; \
	)
	@echo "===== OK: done building doxygen HTML ====="

pdf: $(DOXY_PDF_FILES)
	@echo "===== building doxygen PDF in $(abs_builddir) ====="
	@echo "DOXYINPUT=$(DOXYINPUT)"
	rm -f refman.pdf doxy-pdf*.log
	cd $(abs_builddir) && ( \
	   if doxygen doxy-pdf-config.cfg && ( cd latex && $(MAKE) ); then \
	      $(LN_S) latex/refman.pdf refman.pdf; \
	      exit 0; \
	   else \
	      $(LN_S) latex/refman.log doxy-pdf-make.log; \
	      exit 1; \
	   fi; \
	)
	@echo "===== OK: done building doxygen PDF ====="

clean-local:
	cd $(abs_builddir) && \
	   rm -fr html/ latex/ && \
	   rm -f  refman.pdf doxy-html*.log doxy-pdf*.log
