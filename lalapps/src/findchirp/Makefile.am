## Process this file with automake to produce Makefile.in
# $Id$
if DATAFLOW
DATAFLOWPROGRAMS = lalapps_exc_resp lalapps_trump lalapps_ivana \
                lalapps_icoin lalapps_donald
else
DATAFLOWPROGRAMS = 
endif
bin_PROGRAMS = lalapps_findchirp lalapps_inspinj lalapps_inspinj_s1 \
               lalapps_findchirp_post lalapps_inspawgfile $(DATAFLOWPROGRAMS)
lalapps_findchirp_SOURCES = findchirp.c findchirp.h findchirp_args.c
lalapps_inspinj_SOURCES = inspinj.c
lalapps_inspinj_s1_SOURCES = inspinj_s1.c
lalapps_findchirp_post_SOURCES = findchirp_post.c
lalapps_inspawgfile_SOURCES = inspawgfile.c
if DATAFLOW
lalapps_exc_resp_SOURCES = exc_resp.c
lalapps_trump_SOURCES = trump.c trump.h event_utils.c event_utils.h
lalapps_donald_SOURCES = donald.c donald.h event_utils.c event_utils.h
lalapps_icoin_SOURCES = icoin.c trump.h event_utils.c
lalapps_ivana_SOURCES = ivana.c
endif
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
pkgdata_DATA = BNSMasses.dat
EXTRA_DIST = BNSMasses.dat exc_resp.c trump.c trump.h event_utils.c ivana.c
TESTS_ENVIRONMENT = LALAPPS_DATA_PATH=$(srcdir)
TESTS = lalapps_inspinj_s1
MOSTLYCLEANFILES = injepochs.ilwd injparams.ilwd injlog.txt
