include $(top_srcdir)/gnuscripts/lalapps.am

CC = $(MPICC) -std=gnu99
AM_CPPFLAGS += -I $(srcdir)/..
LIBS += $(MPILIBS)

if LALINFERENCE
if MPI
mpi_progs = lalinference_mcmc
mpi_progs += lalinference_kombine
endif
endif

lalinference_mcmc_SOURCES = \
	LALInferenceMCMC.c \
	LALInferenceMCMCSampler.c

lalinference_kombine_SOURCES = \
	LALInferenceKombine.c \
	LALInferenceKombineSampler.c

bin_PROGRAMS = $(mpi_progs)
