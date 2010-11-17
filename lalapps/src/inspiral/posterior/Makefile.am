SUBDIRS = SPINspiral mpi

LDADD = $(top_builddir)/src/lalapps/liblalapps.la

#WMF Test Programs---remove when time-domain likelihood finished
td_test_progs = time_series_product_test
time_series_product_test_SOURCES = time_series_product_test.c LALInference.c

if LALINSPIRAL
if LALFRAME

if LALMETAIO
metaio_progs = lalapps_inspnest
lalapps_inspnest_SOURCES = lalnest.c LALInspiralMCMC.c LALInspiralMCMCUser.c nest_calc.c nest_calc.h LALInspiralMCMC.h LALInspiralMCMCUser.h
endif

if LALPULSAR
pulsar_progs = lalinference_test lalinference_nest #lalinference_mcmcmpi
lalinference_test_SOURCES = LALInference.c LALInferenceReadData.c LALInferenceTemplate.c InferenceTest.c LALInferenceMCMC.c
lalinference_nest_SOURCES = LALInference.c LALInferenceReadData.c LALInferenceTemplate.c LALInferenceNest.c LALInferenceNestedSampler.c LALInferencePrior.c
#lalinference_mcmcmpi_SOURCES = LALInference.c LALInferenceReadData.c LALInferenceTemplate.c LALInferenceMCMCMPI.c LALInferenceMCMCMPISampler.c LALInferencePrior.c
endif

frame_progs = lalapps_followupMcmc
lalapps_followupMcmc_SOURCES = followupMcmc.c
endif
endif

bin_PROGRAMS = $(metaio_progs) $(frame_progs) $(pulsar_progs)
noinst_PROGRAMS = $(td_test_progs)

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(top_builddir)/src/lalapps
BUILT_SOURCES =
include $(top_srcdir)/misc/vcsID.common

