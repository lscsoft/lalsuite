include $(top_srcdir)/gnuscripts/lalapps.am

MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

SUBDIRS = \
	mpi

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
	lalapps_coherence_test \
	lalapps_nest2pos \
	lalinference_burst_pp_pipe \
	lalinference_pp_pipe \
	lalinference_multi_pipe \
	lalinference_tiger_pipe \
	lalapps_compute_roq_weights \
	lalapps_merge_posteriors
endif

EXTRA_DIST += \
	LALInferenceMCMCSampler.c \
	LALInferenceMCMCSampler.h \
	LALInferenceKombineSampler.c \
	LALInferenceKombineSampler.h
