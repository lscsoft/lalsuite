BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
TESTS =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

SUBDIRS = \
	. \
	mpi \
	$(END_OF_LIST)

# -- C programs -------------

bin_PROGRAMS = \
	lalinference_nest \
	lalinference_injectedlike \
	lalinference_burst \
	lalinference_datadump \
	lalinference_bench \
	lalinference_version \
	$(END_OF_LIST)

lalinference_nest_SOURCES = LALInferenceNest.c
lalinference_injectedlike_SOURCES = LALInferenceInjectedLike.c
lalinference_burst_SOURCES = LALInferenceBurst.c
lalinference_datadump_SOURCES = LALInferenceDataDump.c
lalinference_bench_SOURCES = LALInferenceBench.c
lalinference_version_SOURCES = version.c

TESTS += \
	lalinference_version \
	$(END_OF_LIST)

# -- Bash scripts -----------

dist_bin_SCRIPTS = \
	lalinference_mpi_wrapper \
	$(END_OF_LIST)

# -- Python programs --------

if HAVE_PYTHON

pybin_scripts = \
	cbcBayesBurstPostProc \
	cbcBayesBurstPPAnalysis \
	cbcBayesCombinePosteriors \
	cbcBayesCombinePTMCMCh5s \
	cbcBayesCompPos \
	cbcBayesDIEvidence \
	cbcBayesGraceDBinfo \
	cbcBayesMCMC2pos \
	cbcBayesPlotSpinDisk \
	cbcBayesPosToSimInspiral \
	cbcBayesPosToSimBurst \
	cbcBayesPostProc \
	cbcBayesPPAnalysis \
	cbcBayesThermoInt \
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
	$(END_OF_LIST)

TESTS += \
	$(END_OF_LIST)

endif
