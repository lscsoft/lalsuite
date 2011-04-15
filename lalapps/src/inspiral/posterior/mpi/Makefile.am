LDADD = $(top_builddir)/src/lalapps/liblalapps.la
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(top_builddir)/src/lalapps

CC = $(MPICC)
AM_CFLAGS = -std=gnu99
LIBS += $(MPILIBS)

if MPI
mpi_progs = lalinference_mcmc
lalinference_mcmc_SOURCES = \
	$(srcdir)/../LALInference.c \
	$(srcdir)/../LALInferenceReadData.c \
	$(srcdir)/../LALInferenceTemplate.c \
	$(srcdir)/../LALInferenceMCMC.c \
	$(srcdir)/../LALInferenceMCMCSampler.c \
	$(srcdir)/../LALInferencePrior.c
endif

bin_PROGRAMS = $(mpi_progs)

BUILT_SOURCES =
include $(top_srcdir)/misc/vcsID.common
