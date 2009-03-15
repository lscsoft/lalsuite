if LALFRAME
if LALMETAIO
bin_PROGRAMS = lalapps_inspnest

lalapps_inspnest_SOURCES = lalnest.c LALInspiralMCMC.c LALInspiralMCMCUser.c nest_calc.c nest_calc.h LALInspiralMCMC.h LALInspiralMCMCUser.h
endif
endif

## make sure LALAPPSgitID.h is re-generated at every 'make'
.PHONY: gitID
gitID:
	cd $(top_builddir)/misc && $(MAKE) gitID
BUILT_SOURCES = gitID

LDADD = $(top_builddir)/src/lalapps/liblalapps.la
