python_config.sed:
	@-rm -f python_config.sed
	@echo "s+@PYTHONLIBDIR@+@pkgpythondir@+g" > python_config.sed
	@echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed

laldetchar-hveto: $(srcdir)/laldetchar-hveto.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/laldetchar-hveto.in > laldetchar-hveto
	@chmod +x laldetchar-hveto

bin_SCRIPTS = \
	laldetchar-hveto

EXTRA_DIST = \
	laldetchar-hveto.in

CLEANFILES = \
	laldetchar-hveto \
	python_config.sed
