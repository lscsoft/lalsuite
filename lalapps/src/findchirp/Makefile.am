## Process this file with automake to produce Makefile.in
# $Id$
if LALMETAIO
LALMETAIOPROGRAMS = lalapps_exc_resp lalapps_trump lalapps_ivana \
                lalapps_donald lalapps_snglInspiralReader \
		lalapps_inspinj_find lalapps_gpstogmst lalapps_inspmultiawg
else
LALMETAIOPROGRAMS = 
endif
bin_PROGRAMS = lalapps_inspawgfile lalapps_inspinj $(LALMETAIOPROGRAMS) \
               lalapps_random_bank lalapps_chirplen lalapps_effdist\
               lalapps_BankEfficiency lalapps_BankEfficiencyRealPsd lalapps_InspiralOverlap\
		lalapps_InspiralTools 


lalapps_inspawgfile_SOURCES = inspawgfile.c
lalapps_random_bank_SOURCES = random_bank.c
lalapps_chirplen_SOURCES = chirplen.c
lalapps_effdist_SOURCES = effdist.c
lalapps_inspinj_SOURCES = inspinj.c
lalapps_gpstogmst_SOURCES = gpstogmst.c
lalapps_BankEfficiency_SOURCES = BankEfficiency.c 
lalapps_BankEfficiencyRealPsd_SOURCES = BankEfficiencyRealPsd.c 
lalapps_InspiralOverlap_SOURCES = InspiralOverlap.c 
lalapps_InspiralTools_SOURCES = InspiralTools.c
 
if LALMETAIO
lalapps_exc_resp_SOURCES = exc_resp.c
lalapps_trump_SOURCES = trump.c trump.h event_utils.c event_utils.h
lalapps_donald_SOURCES = donald.c donald.h event_utils.c event_utils.h
lalapps_ivana_SOURCES = ivana.c
lalapps_snglInspiralReader_SOURCES = snglInspiralReader.c event_utils.c event_utils.h
lalapps_inspinj_find_SOURCES = inspinj_find.c event_utils.c event_utils.h
lalapps_inspmultiawg_SOURCES = inspmultiawg.c 
endif
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
pkgdata_DATA = BNSMasses.dat inspsrcs.dat
EXTRA_DIST = BNSMasses.dat inspsrcs.dat exc_resp.c trump.c trump.h event_utils.c ivana.c
MOSTLYCLEANFILES = injepochs.ilwd injparams.ilwd injlog.txt *.xml

