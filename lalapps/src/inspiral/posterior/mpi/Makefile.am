LDADD = $(top_builddir)/src/lalapps/liblalapps.la
AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(top_builddir)/src/lalapps

CC = $(MPICC)
AM_CFLAGS = -std=gnu99
LIBS += $(MPILIBS)

if LALINFERENCE
if MPI
mpi_progs = lalinference_mcmc
endif
endif

lalinference_mcmc_SOURCES = \
	$(srcdir)/../LALInferenceMCMC.c \
	$(srcdir)/../LALInferenceMCMCSampler.c \
	$(srcdir)/../LALInferenceMCMCSampler.h

bin_PROGRAMS = $(mpi_progs)

BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common
