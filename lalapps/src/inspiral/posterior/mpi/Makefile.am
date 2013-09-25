LDADD = $(top_builddir)/src/lalapps/liblalapps.la
AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(top_builddir)/src/lalapps -DPKG_DATA_DIR='"$(pkgdatadir)/"'

CC = $(MPICC)
AM_CFLAGS = -std=gnu99
LIBS += $(MPILIBS)

if LALINFERENCE
if MPI
mpi_progs = lalinference_mcmc
if HAVE_BAMBI
if HAVE_BAMBI_MPI
mpi_progs += lalinference_bambi
AM_CPPFLAGS += -DPARALLEL
AM_CFLAGS += -DPARALLEL
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

BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common
