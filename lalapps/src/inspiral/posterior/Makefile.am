include $(top_srcdir)/gnuscripts/lalapps.am

MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

SUBDIRS = \
	mpi

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
	lalinference_multi_pipe \
	lalinference_tiger_pipe
endif

EXTRA_DIST += \
	LALInferenceMCMCSampler.c \
	LALInferenceMCMCSampler.h \
	LALInferenceKombineSampler.c \
	LALInferenceKombineSampler.h
