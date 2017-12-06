include $(top_srcdir)/gnuscripts/lalapps.am

AM_CPPFLAGS += -I$(top_srcdir)/src/pulsar/Fstatistic

bin_PROGRAMS = lalapps_pulsar_crosscorr lalapps_pulsar_crosscorr_v2

lalapps_pulsar_crosscorr_SOURCES = pulsar_crosscorr.c pulsar_crosscorr.h

lalapps_pulsar_crosscorr_v2_SOURCES = pulsar_crosscorr_v2.c\
        CrossCorrToplist.c CrossCorrToplist.h \
	HeapToplist.c

EXTRA_DIST = testPulsarCrossCorr_v2.sh \
             run_pulsar_crosscorr_v2.py test_Sco_X1.ini \
             combine_crosscorr_toplists.py

TESTS = testPulsarCrossCorr_v2.sh

CLEANFILES = toplist_crosscorr.dat #gsldata.dat
