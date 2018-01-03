BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

SUBDIRS = lalinference

if HAVE_PYTHON

pybin_scripts = \
	bayestar_aggregate_found_injections \
	bayestar_bin_samples \
	bayestar_contour_allsky \
	bayestar_lattice_tmpltbank \
	bayestar_localize_coincs \
	bayestar_localize_lvalert \
	bayestar_plot_allsky \
	bayestar_plot_found_injections \
	bayestar_plot_pileup \
	bayestar_plot_volume \
	bayestar_prune_neighborhood_tmpltbank \
	bayestar_realize_coincs \
	bayestar_sample_model_psd \
	bayestar_sim_to_tmpltbank \
	bayestar_littlehope \
	rapidpe_integrate_extrinsic_likelihood \
	rapidpe_create_event_dag \
	rapidpe_compute_intrinsic_fisher \
	rapidpe_compute_intrinsic_grid \
	rapidpe_calculate_overlap \
	imrtgr_imr_consistency_test \
	lalinference_pipe \
	$(END_OF_LIST)

endif
