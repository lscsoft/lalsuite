python_config.sed:
	@-rm -f python_config.sed
	@echo "s+@PYTHONLIBDIR@+@pkgpythondir@+g" > python_config.sed
	@echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed

laldetchar-hveto: $(srcdir)/laldetchar-hveto.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/laldetchar-hveto.in > laldetchar-hveto
	@chmod +x laldetchar-hveto

laldetchar-hveto-summary-page: $(srcdir)/laldetchar-hveto-summary-page.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/laldetchar-hveto-summary-page.in > laldetchar-hveto-summary-page
	@chmod +x laldetchar-hveto-summary-page

laldetchar-segments-from-data: $(srcdir)/laldetchar-segments-from-data.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/laldetchar-segments-from-data.in > laldetchar-segments-from-data
	@chmod +x laldetchar-segments-from-data

bin_SCRIPTS = \
	laldetchar-hveto \
	laldetchar-hveto-summary-page \
	laldetchar-segments-from-data

EXTRA_DIST = \
	laldetchar-hveto.in \
	laldetchar-hveto-summary-page.in \
	laldetchar-segments-from-data.in

CLEANFILES = \
	laldetchar-hveto \
	laldetchar-hveto-summary-page \
	laldetchar-segments-from-data \
	python_config.sed
