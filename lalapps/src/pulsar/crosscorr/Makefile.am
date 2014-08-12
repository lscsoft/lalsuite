BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/lalapps.am

AM_CPPFLAGS = -I$(top_srcdir)/src/pulsar/FDS_isolated @AM_CPPFLAGS@

CFSsrc = $(top_srcdir)/src/pulsar/FDS_isolated

bin_PROGRAMS = lalapps_pulsar_crosscorr lalapps_pulsar_crosscorr_v2

lalapps_pulsar_crosscorr_SOURCES = pulsar_crosscorr.c pulsar_crosscorr.h

lalapps_pulsar_crosscorr_v2_SOURCES = pulsar_crosscorr_v2.c\
        CrossCorrToplist.c CrossCorrToplist.h \
	../FDS_isolated/HeapToplist.c ../FDS_isolated/HeapToplist.h

EXTRA_DIST = testPulsarCrossCorr_v2.sh \
             run_pulsar_crosscorr_v2.py test_Sco_X1.ini

TESTS = testPulsarCrossCorr_v2.sh

TESTS_ENVIRONMENT = \
	LAL_DEBUG_LEVEL=msglvl1,memdbg; export LAL_DEBUG_LEVEL; \
	LAL_DATA_PATH="$(LAL_DATA_PATH)"; export LAL_DATA_PATH;

CLEANFILES = toplist_crosscorr.dat SFTpairlist.dat SFTlist.dat