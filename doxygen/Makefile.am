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
                          $(addprefix $(abs_top_srcdir)/,$(DOXYINPUT))

override TEXINPUTS = .:$(abs_srcdir):
override pool_size = 4000000

export doxy_inputdirs TEXINPUTS pool_size

html: doxygen-html.cfg
	@echo "***** building doxygen HTML in $(abs_builddir) *****"
	@echo "DOXYINPUT=$(DOXYINPUT)"
	doxygen $(abs_builddir)/doxygen-html.cfg
	@echo "***** done building doxygen HTML *****"

pdf: doxygen-pdf.cfg
	@echo "***** building doxygen PDF in $(abs_builddir) *****"
	@echo "DOXYINPUT=$(DOXYINPUT)"
	doxygen $(abs_builddir)/doxygen-pdf.cfg
	cd $(abs_builddir)/latex && $(MAKE)
	cd $(abs_builddir) && rm -f refman.pdf && $(LN_S) latex/refman.pdf .
	@echo "***** done building doxygen PDF *****"

clean-local:
	cd $(abs_builddir) && \
	rm -fr html/ latex/ && \
	rm -f  refman.pdf doxygen-html.log doxygen-pdf.log
