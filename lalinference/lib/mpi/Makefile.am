MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalsuite_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

LDADD += ../liblalinference.la

CC = $(MPICC) -std=gnu99
LIBS += $(MPILIBS)

if MPI
mpi_progs = lalinference_mcmc
mpi_progs += lalinference_kombine
endif

lalinference_mcmc_SOURCES = \
	LALInferenceMCMC.c \
	LALInferenceMCMCSampler.c

lalinference_kombine_SOURCES = \
	LALInferenceKombine.c \
	LALInferenceKombineSampler.c

noinst_HEADERS = \
	LALInferenceKombineSampler.h \
	LALInferenceMCMCSampler.h \
	$(END_OF_LIST)

bin_PROGRAMS = $(mpi_progs)
