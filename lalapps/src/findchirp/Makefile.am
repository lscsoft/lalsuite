## Process this file with automake to produce Makefile.in
# $Id$
if LALMETAIO
LALMETAIOPROGRAMS = lalapps_exc_resp lalapps_ivana lalapps_inspmultiawg \
                    lalapps_chirplen 
else
LALMETAIOPROGRAMS = 
endif
bin_PROGRAMS = lalapps_inspawgfile $(LALMETAIOPROGRAMS) \
               lalapps_random_bank  lalapps_effdist\
               lalapps_InspiralOverlap lalapps_inspiralDistance lalapps_DistanceVsMass \
               lalapps_BankNumberOfTemplates 


lalapps_DistanceVsMass_SOURCES 		= DistanceVsMass.c
lalapps_inspiralDistance_SOURCES 	= inspiralDistance.c
lalapps_inspawgfile_SOURCES 		= inspawgfile.c
lalapps_random_bank_SOURCES 		= random_bank.c
lalapps_effdist_SOURCES 		= effdist.c
lalapps_BankNumberOfTemplates_SOURCES	= BankNumberOfTemplates.c
lalapps_InspiralOverlap_SOURCES 	= InspiralOverlap.c 

if LALMETAIO
lalapps_exc_resp_SOURCES 		= exc_resp.c
lalapps_ivana_SOURCES 			= ivana.c
lalapps_inspmultiawg_SOURCES 		= inspmultiawg.c 
lalapps_chirplen_SOURCES 		= chirplen.c
lalapps_BankEfficiency_SOURCES		= BankEfficiency.c
endif
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

EXTRA_DIST = 	exc_resp.c ivana.c inspmultiawg.c chirplen.c \
		BE_create_script.py BEAscii2xml.c BE_process_xml.py

CLEANFILES = 	lalapps_BankEfficiency  \
		lalapps_BankNumberOfTemplates  lalapps_InspiralOverlap

