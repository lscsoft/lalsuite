BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/lalapps.am

if LALINSPIRAL
LALINSPIRALPROGRAMS = \
	lalapps_InspiralOverlap \
	lalapps_inspiralDistance \
	lalapps_DistanceVsMass \
	lalapps_random_bank \
	lalapps_effdist \
	lalapps_inspawgfile
if LALMETAIO
LALMETAIOPROGRAMS = \
	lalapps_exc_resp \
	lalapps_ivana \
	lalapps_chirplen \
	lalapps_BankEfficiency
endif
endif
bin_PROGRAMS = $(LALINSPIRALPROGRAMS) $(LALMETAIOPROGRAMS)

if LALINSPIRAL
lalapps_DistanceVsMass_SOURCES = DistanceVsMass.c
lalapps_inspiralDistance_SOURCES = inspiralDistance.c
lalapps_inspawgfile_SOURCES = inspawgfile.c
lalapps_random_bank_SOURCES = random_bank.c
lalapps_effdist_SOURCES = effdist.c
lalapps_InspiralOverlap_SOURCES = InspiralOverlap.c
if LALMETAIO
lalapps_exc_resp_SOURCES = exc_resp.c
lalapps_ivana_SOURCES = ivana.c
lalapps_chirplen_SOURCES = chirplen.c
lalapps_BankEfficiency_SOURCES = BankEfficiency.c BankEfficiency.h
endif
endif

EXTRA_DIST = \
	exc_resp.c \
	ivana.c \
	chirplen.c \
	TestRandomBank.c \
	BankEfficiency.c \
	BankEfficiency.h \
	BE_create_script.py \
	TestRandomBank.h

CLEANFILES = \
	lalapps_BankEfficiency  \
	lalapps_InspiralOverlap
