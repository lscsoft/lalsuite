BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

SUBDIRS = lalinference

if HAVE_PYTHON

pybin_scripts = \
	rapidpe_integrate_extrinsic_likelihood \
	rapidpe_create_event_dag \
	rapidpe_compute_intrinsic_fisher \
	rapidpe_compute_intrinsic_grid \
	rapidpe_calculate_overlap \
	imrtgr_imr_consistency_test \
	lalinference_burst_pp_pipe \
	lalinference_coherence_test \
	lalinference_compute_roq_weights \
	lalinference_merge_posteriors \
	lalinference_multi_pipe \
	lalinference_nest2pos \
	lalinference_pipe \
	lalinference_pp_pipe \
	lalinference_review_test \
	lalinference_tiger_pipe \
	cbcBayesBurstPostProc \
	cbcBayesCombinePosteriors \
	cbcBayesCompPos \
	cbcBayesDIEvidence \
	cbcBayesPlotSpinDisk \
	cbcBayesPosToSimInspiral \
	cbcBayesPPAnalysis \
	cbcBayesBurstPPAnalysis \
	cbcBayesCombinePTMCMCh5s \
	cbcBayesMCMC2pos \
	cbcBayesPosToSimBurst \
	cbcBayesPostProc \
	cbcBayesThermoInt \
	$(END_OF_LIST)

endif
