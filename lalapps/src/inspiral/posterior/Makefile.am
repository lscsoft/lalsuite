include $(top_srcdir)/gnuscripts/lalapps.am

MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

SUBDIRS = \
	mpi

if LALPULSAR
if LALINFERENCE

pulsar_progs = lalinference_nest \
	lalinference_bench \
	lalinference_datadump \
	lalinference_burst
lalinference_nest_SOURCES = LALInferenceNest.c
lalinference_burst_SOURCES = LALInferenceBurst.c
lalinference_datadump_SOURCES = LALInferenceDataDump.c
lalinference_bench_SOURCES = LALInferenceBench.c
if HAVE_BAMBI
if !HAVE_BAMBI_MPI
pulsar_progs += lalinference_bambi
lalinference_bambi_SOURCES = LALInferenceBAMBI.c
lalinference_bambi_LDADD = ${BAMBI_LIBS}
endif # bambi_mpi
endif # bambi

endif # lalinference
endif # lalpulsar

#bayesline_SOURCES = \
#	BayesLine.c

dist_bin_SCRIPTS = \
	lalinference_mpi_wrapper

bin_PROGRAMS = \
	$(metaio_progs) \
	$(frame_progs) \
	$(pulsar_progs) \
	$(mpi_progs)

noinst_PROGRAMS = \
	$(td_test_progs)

if HAVE_PYTHON
pybin_scripts = \
	lalapps_combine_evidence \
	lalapps_coherence_test \
	lalapps_nest_pipe \
	lalapps_vo_nest2pos \
	lalinference_pipe \
	lalapps_nest2pos \
	lalapps_merge_nested_sampling_runs \
	lalinference_burst_pp_pipe \
	lalinference_pp_pipe \
	lalinference_multi_pipe \
	lalinference_tiger_pipe \
	lalapps_compute_roq_weights \
	lalapps_merge_posteriors \
	lalapps_test_run_lalinference_pipe_example
pkgpython_PYTHON = \
	nest_utils.py \
	combine_evidence.py \
	nest2pos.py \
	lalinference_pipe_utils.py
endif

EXTRA_DIST += \
	LALInferenceMCMCSampler.c \
	LALInferenceMCMCSampler.h \
	LALInferenceKombineSampler.c \
	LALInferenceKombineSampler.h
