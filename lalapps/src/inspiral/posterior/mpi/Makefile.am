include $(top_srcdir)/gnuscripts/lalapps.am

CC = $(MPICC) -std=gnu99
AM_CPPFLAGS += -I $(srcdir)/..
LIBS += $(MPILIBS)

if LALINFERENCE
if MPI
mpi_progs = lalinference_mcmc
mpi_progs += lalinference_kombine
if HAVE_BAMBI
if HAVE_BAMBI_MPI
mpi_progs += lalinference_bambi
AM_CPPFLAGS += -DPARALLEL
lalinference_bambi_SOURCES = LALInferenceBAMBI.c
lalinference_bambi_LDADD = ${BAMBI_LIBS}
lalinference_bambi_LDFLAGS = $(shell ${MPICXX} --showme:link) $(shell ${MPIFC} --showme:link)
endif
endif
endif
endif

lalinference_mcmc_SOURCES = \
	LALInferenceMCMC.c \
	LALInferenceMCMCSampler.c

lalinference_kombine_SOURCES = \
	LALInferenceKombine.c \
	LALInferenceKombineSampler.c

bin_PROGRAMS = $(mpi_progs)
