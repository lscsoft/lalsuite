EXTRA_DIST = \
	DoxygenLayout-html.xml \
	DoxygenLayout-pdf.xml \
	doxy-html-aliases.cfg \
	doxy-pdf-aliases.cfg \
	doxy_html_eqns.sty \
	doxy_pdf_header.tex \
	doxygen-html.cfg.in \
	doxygen-pdf.cfg.in

SUBDIRS = src

all: html

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
      lsv="$$lsv-"`$(GIT) log -1 --format="%h"`; \
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

html: doxygen-html.cfg
	@echo "===== building doxygen HTML in $(abs_builddir) ======"
	@echo "DOXYINPUT=$(DOXYINPUT)"
	rm -f doxygen-html*.log
	cd $(abs_builddir) && ( \
	   if doxygen doxygen-html.cfg; then \
	      if test -f html/_formulas.log; then \
	         $(LN_S) html/_formulas.log doxygen-html-eqns.log; \
	      fi; \
	      exit 0; \
	   fi; \
	)
	@echo "===== OK: done building doxygen HTML ====="

pdf: doxygen-pdf.cfg
	@echo "===== building doxygen PDF in $(abs_builddir) ====="
	@echo "DOXYINPUT=$(DOXYINPUT)"
	rm -f refman.pdf doxygen-pdf*.log
	cd $(abs_builddir) && ( \
	   if doxygen doxygen-pdf.cfg && ( cd latex && $(MAKE) ); then \
	      $(LN_S) latex/refman.pdf refman.pdf; \
	      exit 0; \
	   else \
	      $(LN_S) latex/refman.log doxygen-pdf-make.log; \
	      exit 1; \
	   fi; \
	)
	@echo "===== OK: done building doxygen PDF ====="

clean-local:
	cd $(abs_builddir) && \
	   rm -fr html/ latex/ && \
	   rm -f  refman.pdf doxygen-html*.log doxygen-pdf*.log
