if LALFRAME
if LALMETAIO
if LALINSPIRAL
bin_PROGRAMS = lalapps_spinspiral
lalapps_spinspiral_SOURCES = \
	SPINspiral_main.c \
	SPINspiral_mcmc.c \
	SPINspiral_templates.c \
	SPINspiral_routines.c \
	SPINspiral_data.c \
	SPINspiral_parameters.c \
	SPINspiral_signal.c \
	SPINspiral_lal.c \
	SPINspiral_3rdparty.c \
	SPINspiral.h remez.h
endif
endif
endif

LDADD = $(top_builddir)/src/lalapps/liblalapps.la

clean-local :
	rm -f mcmc.output.??????.?? *-signal*.dat
