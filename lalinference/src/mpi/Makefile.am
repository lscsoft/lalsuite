include $(top_srcdir)/gnuscripts/lalsuite_header_links.am

CC = $(MPICC) -std=gnu99
AM_CPPFLAGS += -I $(srcdir)/..
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

bin_PROGRAMS = $(mpi_progs)
