## make sure LALAPPSgitID.h is re-generated at every 'make'
.PHONY: gitID
gitID:
	cd $(top_builddir)/misc && $(MAKE) gitID
BUILT_SOURCES = gitID

SUBDIRS = SPINspiral

LDADD = $(top_builddir)/src/lalapps/liblalapps.la


if LALFRAME

if LALMETAIO
metaio_progs = lalapps_inspnest
lalapps_inspnest_SOURCES = lalnest.c LALInspiralMCMC.c LALInspiralMCMCUser.c nest_calc.c nest_calc.h LALInspiralMCMC.h LALInspiralMCMCUser.h
endif

frame_progs = lalapps_followupMcmc
lalapps_followupMcmc_SOURCES = followupMcmc.c

endif

bin_PROGRAMS = $(metaio_progs) $(frame_progs)
