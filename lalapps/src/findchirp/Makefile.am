## Process this file with automake to produce Makefile.in
# $Id$
if LALMETAIO
LALMETAIOPROGRAMS = lalapps_exc_resp lalapps_ivana lalapps_inspmultiawg
else
LALMETAIOPROGRAMS = 
endif
bin_PROGRAMS = lalapps_inspawgfile $(LALMETAIOPROGRAMS) \
               lalapps_random_bank lalapps_chirplen lalapps_effdist\
               lalapps_BankEfficiency lalapps_BankEfficiencyRealPsd \
               lalapps_InspiralOverlap lalapps_InspiralTools  


lalapps_inspawgfile_SOURCES = inspawgfile.c
lalapps_random_bank_SOURCES = random_bank.c
lalapps_chirplen_SOURCES = chirplen.c
lalapps_effdist_SOURCES = effdist.c
lalapps_BankEfficiency_SOURCES = BankEfficiency.c 
lalapps_BankEfficiencyRealPsd_SOURCES = BankEfficiencyRealPsd.c 
lalapps_InspiralOverlap_SOURCES = InspiralOverlap.c 
lalapps_InspiralTools_SOURCES = InspiralTools.c
 
if LALMETAIO
lalapps_exc_resp_SOURCES = exc_resp.c
lalapps_ivana_SOURCES = ivana.c
lalapps_inspmultiawg_SOURCES = inspmultiawg.c 
endif
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
EXTRA_DIST = exc_resp.c ivana.c

