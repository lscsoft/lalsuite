## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

SUBDIRS =  MakeSFTs Injections FDS_isolated hough GCT FDS_binary SideBand templateBanks crosscorr TimingTests TDS_isolated fscan TwoSpect CombSearch SFTReferenceLibrary

if LALFRAME
# this only requires LALFrame because of usage of LALFrameL.h
LALFRAMEPROGS = lalapps_make_sfts
else
LALFRAMEPROGS =
endif

bin_PROGRAMS = $(LALFRAMEPROGS) lalapps_Calibrate60sSFTs lalapps_Calibrate1800SFTs lalapps_FreqAverager lalapps_psd \
	lalapps_ComputePSD  lalapps_FreqAverager_v2 lalapps_GetSFTVelocities lalapps_create_solar_system_ephemeris \
        lalapps_make_job_input_data lalapps_PrintDetectorState lalapps_create_time_correction_ephemeris

if LALFRAME
lalapps_make_sfts_SOURCES = make_sfts.c
endif

lalapps_Calibrate60sSFTs_SOURCES = Calibrate60sSFTs.c Calibrate60sSFTs.h
lalapps_Calibrate1800SFTs_SOURCES = Calibrate1800SFTs.c Calibrate1800SFTs.h
lalapps_FreqAverager_SOURCES = FreqAverager.c FreqAverager.h
lalapps_psd_SOURCES = psd.c psd.h
lalapps_ComputePSD_SOURCES = ComputePSD.c
lalapps_FreqAverager_v2_SOURCES =  FreqAverager_v2.c
lalapps_GetSFTVelocities_SOURCES = GetSFTVelocities.c
lalapps_create_solar_system_ephemeris_SOURCES = create_solar_system_ephemeris.c
lalapps_make_job_input_data_SOURCES = make_job_input_data.c
lalapps_PrintDetectorState_SOURCES = PrintDetectorState.c
lalapps_create_time_correction_ephemeris_SOURCES = create_time_correction_ephemeris.c create_time_correction_ephemeris.h

EXTRA_DIST = \
make_locked_data_list \
make_cluster_data_list \
make_job_input_data.c \
make_sfts.c \
Calibrate60sSFTs.c \
Calibrate60sSFTs.h \
Calibrate1800SFTs.c \
Calibrate1800SFTs.h \
FreqAverager.c \
FreqAverager.h \
psd.c \
psd.h \
README \
SFT.00000 \
SFT.00001 \
SFT.00002 \
SFT.00003 \
SFT.00004 \
SFT.00005 \
SFT.00006 \
SFT.00007 \
SFT.00008 \
SFT.00009 \
S1-H1segments.txt \
S1-H2segments.txt \
S1-L1segments.txt \
S2-H1segments.txt \
S2-H2segments.txt \
S2-L1segments.txt \
testComputePSD.sh \
psd_ref.dat


TESTS = testComputePSD.sh

CLEANFILES = psd1.dat

TESTS_ENVIRONMENT = 1>/dev/null\
LALSUITE_TOP_SRCDIR=$(LALSUITE_TOP_SRCDIR) \
LAL_DATA_PATH=$(LAL_DATA_PATH):.:$${LALPULSAR_PREFIX:+"$${LALPULSAR_PREFIX}/share/lalpulsar"}:$${LALSUITE_TOP_SRCDIR:+"$${LALSUITE_TOP_SRCDIR}/lalpulsar/test/"}

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/misc/vcsID.common
