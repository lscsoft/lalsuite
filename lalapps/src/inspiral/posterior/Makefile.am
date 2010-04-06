SUBDIRS = SPINspiral

LDADD = $(top_builddir)/src/lalapps/liblalapps.la

if LALINSPIRAL
if LALFRAME

if LALMETAIO
metaio_progs = lalapps_inspnest lalinference_test
lalapps_inspnest_SOURCES = lalnest.c LALInspiralMCMC.c LALInspiralMCMCUser.c nest_calc.c nest_calc.h LALInspiralMCMC.h LALInspiralMCMCUser.h
endif
lalinference_test_SOURCES = LALInference.c LALInferenceReadData.c LALInferenceTemplate.c InferenceTest.c

frame_progs = lalapps_followupMcmc
lalapps_followupMcmc_SOURCES = followupMcmc.c
endif
endif

bin_PROGRAMS = $(metaio_progs) $(frame_progs)
