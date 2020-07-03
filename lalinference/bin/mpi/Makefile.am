MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalsuite_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

LDADD += ../../lib/liblalinference.la

CC = $(MPICC) -std=gnu99
LIBS += $(MPILIBS)

if MPI

bin_PROGRAMS = \
	lalinference_mcmc \
	lalinference_kombine \
	$(END_OF_LIST)

lalinference_mcmc_SOURCES = \
	LALInferenceMCMC.c \
	LALInferenceMCMCSampler.c \
	$(END_OF_LIST)

lalinference_kombine_SOURCES = \
	LALInferenceKombine.c \
	LALInferenceKombineSampler.c \
	$(END_OF_LIST)

noinst_HEADERS = \
	LALInferenceKombineSampler.h \
	LALInferenceMCMCSampler.h \
	$(END_OF_LIST)

endif
