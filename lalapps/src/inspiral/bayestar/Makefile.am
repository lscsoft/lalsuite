## Process this file with automake to produce Makefile.in

if HAVE_PYTHON
bin_SCRIPTS = \
	bayestar_aggregate_found_injections \
	bayestar_lattice_tmpltbank \
	bayestar_localize_coincs \
	bayestar_localize_lvalert \
	bayestar_plot_allsky \
	bayestar_plot_found_injections \
	bayestar_plot_pileup \
	bayestar_prune_neighborhood_tmpltbank \
	bayestar_realize_coincs \
	bayestar_sim_to_tmpltbank \
	bayestar_littlehope
endif

python_config.sed: Makefile
	@-rm -f python_config.sed
	@echo "s+@PYTHONLIBDIR@+@pkgpythondir@+g" > python_config.sed
	@echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed

bayestar_aggregate_found_injections: bayestar_aggregate_found_injections.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $< > $@

bayestar_lattice_tmpltbank: bayestar_lattice_tmpltbank.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $< > $@

bayestar_localize_coincs: bayestar_localize_coincs.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $< > $@

bayestar_localize_lvalert: bayestar_localize_lvalert.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $< > $@

bayestar_plot_allsky: bayestar_plot_allsky.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $< > $@

bayestar_plot_found_injections: bayestar_plot_found_injections.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $< > $@

bayestar_plot_pileup: bayestar_plot_pileup.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $< > $@

bayestar_prune_neighborhood_tmpltbank: bayestar_prune_neighborhood_tmpltbank.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $< > $@

bayestar_realize_coincs: bayestar_realize_coincs.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $< > $@

bayestar_sim_to_tmpltbank: bayestar_sim_to_tmpltbank.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $< > $@

bayestar_littlehope: bayestar_littlehope.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $< > $@

EXTRA_DIST = \
	bayestar_aggregate_found_injections.in \
	bayestar_lattice_tmpltbank.in \
	bayestar_localize_coincs.in \
	bayestar_localize_lvalert.in \
	bayestar_plot_allsky.in \
	bayestar_plot_found_injections.in \
	bayestar_plot_pileup.in \
	bayestar_prune_neighborhood_tmpltbank.in \
	bayestar_realize_coincs.in \
	bayestar_sim_to_tmpltbank.in \
	bayestar_littlehope.in

CLEANFILES = \
	$(bin_SCRIPTS) \
	python_config.sed
