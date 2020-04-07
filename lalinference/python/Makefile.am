BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

SUBDIRS = lalinference

if HAVE_PYTHON

pybin_scripts = \
	imrtgr_imr_consistency_test \
	lalinference_burst_pp_pipe \
	lalinference_coherence_test \
	lalinference_compute_roq_weights \
	lalinference_cpnest \
	lalinference_evolve_spins_and_append_samples \
	lalinference_merge_posteriors \
	lalinference_multi_pipe \
	lalinference_nest2pos \
	lalinference_pipe \
	lalinference_pp_pipe \
	lalinference_review_test \
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
	cbcBayesGraceDBinfo \
	cbcBayesPosToSimBurst \
	cbcBayesPostProc \
	cbcBayesThermoInt \
	$(END_OF_LIST)

endif
