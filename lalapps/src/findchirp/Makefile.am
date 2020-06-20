BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

LALINSPIRALPROGRAMS = \
	lalapps_inspiralDistance \
	lalapps_DistanceVsMass \
	lalapps_random_bank \
	lalapps_effdist \
	lalapps_inspawgfile
if LALMETAIO
LALMETAIOPROGRAMS = \
	lalapps_exc_resp \
	lalapps_chirplen
endif
bin_PROGRAMS = $(LALINSPIRALPROGRAMS) $(LALMETAIOPROGRAMS)

lalapps_DistanceVsMass_SOURCES = DistanceVsMass.c
lalapps_inspiralDistance_SOURCES = inspiralDistance.c
lalapps_inspawgfile_SOURCES = inspawgfile.c
lalapps_random_bank_SOURCES = random_bank.c
lalapps_effdist_SOURCES = effdist.c
if LALMETAIO
lalapps_exc_resp_SOURCES = exc_resp.c
lalapps_chirplen_SOURCES = chirplen.c
endif

EXTRA_DIST = \
	chirplen.c \
	exc_resp.c \
	$(END_OF_LIST)
