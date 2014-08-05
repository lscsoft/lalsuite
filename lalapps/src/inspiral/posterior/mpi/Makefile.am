BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/lalapps.am

CC = $(MPICC) -std=gnu99
AM_CPPFLAGS = @AM_CPPFLAGS@
LIBS += $(MPILIBS)

if LALINFERENCE
if MPI
mpi_progs = lalinference_mcmc
if HAVE_BAMBI
if HAVE_BAMBI_MPI
mpi_progs += lalinference_bambi
AM_CPPFLAGS += -DPARALLEL
lalinference_bambi_SOURCES = $(srcdir)/../LALInferenceBAMBI.c
lalinference_bambi_LDADD = ${BAMBI_LIBS}
lalinference_bambi_LDFLAGS = $(shell ${MPICXX} --showme:link) $(shell ${MPIFC} --showme:link)
endif
endif
endif
endif

lalinference_mcmc_SOURCES = \
	$(srcdir)/../LALInferenceMCMC.c \
	$(srcdir)/../LALInferenceMCMCSampler.c \
	$(srcdir)/../LALInferenceMCMCSampler.h

bin_PROGRAMS = $(mpi_progs)
