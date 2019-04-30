include $(top_srcdir)/gnuscripts/lalsuite_header_links.am

BUILT_SOURCES = liblalinference

.PHONY: liblalinference
liblalinference:
	@cd $(top_builddir)/src && $(MAKE) liblalinference.la

#LDADD += $(top_builddir)/src/liblalinference.la

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
