all: html

.PHONY: html pdf

export TEXINPUTS = .:$(abs_srcdir):
export pool_size = 4000000

html:
	@echo "$(abs_builddir)/Makefile: Running doxygen to produce HTML ..."
	doxygen $(abs_builddir)/doxygen-html.cfg
	@echo "$(abs_builddir)/Makefile: ... done."

pdf:
	@echo "$(abs_builddir)/Makefile: Running doxygen to produce PDF ..."
	doxygen $(abs_builddir)/doxygen-pdf.cfg
	@echo "$(abs_builddir)/Makefile: Compiling the PDF manual 'refman.pdf' ..."
	cd $(abs_builddir)/latex && $(MAKE)
	cd $(abs_builddir) && rm -f refman.pdf && $(LN_S) latex/refman.pdf .
	@echo "$(abs_builddir)/Makefile: ... done."

clean-local:
	cd $(abs_builddir) && rm -rf html/ latex/ refman.pdf
