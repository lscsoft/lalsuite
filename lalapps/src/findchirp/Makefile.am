## Process this file with automake to produce Makefile.in
# $Id$
if DATAFLOW
DATAFLOWPROGRAMS = lalapps_exc_resp lalapps_trump lalapps_ivana \
                lalapps_donald lalapps_snglInspiralReader \
                lalapps_inspinj lalapps_inspinj_s1 lalapps_inca
else
DATAFLOWPROGRAMS = 
endif
bin_PROGRAMS = lalapps_findchirp_post lalapps_inspawgfile $(DATAFLOWPROGRAMS) \
               lalapps_random_bank lalapps_chirplen lalapps_effdist
lalapps_findchirp_post_SOURCES = findchirp_post.c
lalapps_inspawgfile_SOURCES = inspawgfile.c
lalapps_random_bank_SOURCES = random_bank.c
lalapps_chirplen_SOURCES = chirplen.c
lalapps_effdist_SOURCES = effdist.c
if DATAFLOW
lalapps_inspinj_SOURCES = inspinj.c
lalapps_inca_SOURCES = inca.c event_utils.c event_utils.h
lalapps_inspinj_s1_SOURCES = inspinj_s1.c
lalapps_exc_resp_SOURCES = exc_resp.c
lalapps_trump_SOURCES = trump.c trump.h event_utils.c event_utils.h
lalapps_donald_SOURCES = donald.c donald.h event_utils.c event_utils.h
lalapps_ivana_SOURCES = ivana.c
lalapps_snglInspiralReader_SOURCES = snglInspiralReader.c event_utils.c event_utils.h
endif
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
pkgdata_DATA = BNSMasses.dat inspsrcs.dat
EXTRA_DIST = BNSMasses.dat inspsrcs.dat exc_resp.c trump.c trump.h event_utils.c ivana.c
MOSTLYCLEANFILES = injepochs.ilwd injparams.ilwd injlog.txt
