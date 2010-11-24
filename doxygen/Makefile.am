all: html

.PHONY: getfigs html pdf

export abs_builddir abs_srcdir abs_lalsuite_srcdir MKDIR_P LN_S

export TEXINPUTS = .:$(abs_srcdir):
export pool_size = 4000000

getfigs:
	@echo "$(abs_builddir)/Makefile: Collecting figures from LALSuite ..."
	$(abs_srcdir)/collect_figures.sh
	@echo "$(abs_builddir)/Makefile: ... done."

html: getfigs
	@echo "$(abs_builddir)/Makefile: Running doxygen to produce HTML ..."
	doxygen $(abs_builddir)/doxygen-html.cfg
	@echo "$(abs_builddir)/Makefile: ... done."

pdf: getfigs
	@echo "$(abs_builddir)/Makefile: Running doxygen to produce PDF ..."
	doxygen $(abs_builddir)/doxygen-pdf.cfg
	@echo "$(abs_builddir)/Makefile: Compiling the PDF manual 'refman.pdf' ..."
	cd $(abs_builddir)/latex && $(MAKE)
	cd $(abs_builddir) && rm -f refman.pdf && $(LN_S) latex/refman.pdf .
	@echo "$(abs_builddir)/Makefile: ... done."

clean-local:
	cd $(abs_builddir) && rm -rf html latex figures refman.pdf
