BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalapps.am

SUBDIRS = \
	SPINspiral \
	mpi

if LALINSPIRAL
if LALFRAME

if LALPULSAR
if LALINFERENCE
pulsar_progs = lalinference_nest ChiSquare_test #lalinference_mcmcmpi
ChiSquare_test_SOURCES = ChiSquare_test.c
lalinference_nest_SOURCES = LALInferenceNest.c
if HAVE_BAMBI
if !HAVE_BAMBI_MPI
pulsar_progs += lalinference_bambi
lalinference_bambi_SOURCES = LALInferenceBAMBI.c
lalinference_bambi_LDADD = ${BAMBI_LIBS}
endif # bambi_mpi
endif # bambi

endif # lalinference
endif # lalpulsar

frame_progs = \
	lalapps_followupMcmc
lalapps_followupMcmc_SOURCES = \
	followupMcmc.c

endif # lalframe
endif # lalinspiral

#bayesline_SOURCES = \
#	BayesLine.c

bin_PROGRAMS = \
	$(metaio_progs) \
	$(frame_progs) \
	$(pulsar_progs) \
	$(mpi_progs)
#	bayesline

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
	lalinference_pp_pipe \
	lalapps_compute_roq_weights
pkgpython_PYTHON = \
	nest_utils.py \
	combine_evidence.py \
	nest2pos.py \
	lalinference_pipe_utils.py
endif

EXTRA_DIST += \
	LALInferenceMCMCSampler.c \
	LALInferenceMCMCSampler.h \
	LALInferenceEnsembleSampler.c \
	LALInferenceEnsembleSampler.h
