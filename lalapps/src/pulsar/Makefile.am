## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(top_builddir)/src/lalapps -DPKG_DATA_DIR='"$(pkgdatadir)/"' @AM_CPPFLAGS@
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

SUBDIRS = \
	MakeSFTs \
	Injections \
	FDS_isolated \
	hough \
	GCT \
	FDS_binary \
	SideBand \
	templateBanks \
	crosscorr \
	TimingTests \
	TDS_isolated \
	fscan \
	TwoSpect \
	CombSearch \
	SFTReferenceLibrary \
	CreateEphemeris

bin_PROGRAMS = \
	lalapps_Calibrate60sSFTs \
	lalapps_Calibrate1800SFTs \
	lalapps_FreqAverager \
	lalapps_psd \
	lalapps_ComputePSD \
	lalapps_FreqAverager_v2 \
	lalapps_GetSFTVelocities \
	lalapps_make_job_input_data \
	lalapps_PrintDetectorState \
	lalapps_ComputeAntennaPattern

if LALFRAME
# this only requires LALFrame because of usage of LALFrameL.h
bin_PROGRAMS += lalapps_make_sfts
lalapps_make_sfts_SOURCES = make_sfts.c
endif

lalapps_Calibrate60sSFTs_SOURCES = Calibrate60sSFTs.c
lalapps_Calibrate1800SFTs_SOURCES = Calibrate1800SFTs.c
lalapps_FreqAverager_SOURCES = FreqAverager.c
lalapps_psd_SOURCES = psd.c psd.h
lalapps_ComputePSD_SOURCES = ComputePSD.c
lalapps_FreqAverager_v2_SOURCES =  FreqAverager_v2.c
lalapps_GetSFTVelocities_SOURCES = GetSFTVelocities.c
lalapps_make_job_input_data_SOURCES = make_job_input_data.c
lalapps_PrintDetectorState_SOURCES = PrintDetectorState.c
lalapps_ComputeAntennaPattern_SOURCES = ComputeAntennaPattern.c

EXTRA_DIST = \
	make_locked_data_list \
	make_cluster_data_list \
	make_job_input_data.c \
	make_sfts.c \
	Calibrate60sSFTs.c \
	Calibrate1800SFTs.c \
	FreqAverager.c \
	psd.c \
	psd.h \
	README \
	S1-H1segments.txt \
	S1-H2segments.txt \
	S1-L1segments.txt \
	S2-H1segments.txt \
	S2-H2segments.txt \
	S2-L1segments.txt \
	testComputePSD.sh \
	testComputeAntennaPattern.sh

TESTS = \
	testComputePSD.sh
#	testComputeAntennaPattern.sh

TESTS_ENVIRONMENT = \
	LAL_DEBUG_LEVEL=memdbg; export LAL_DEBUG_LEVEL; \
	LAL_DATA_PATH="$(LAL_DATA_PATH)"; export LAL_DATA_PATH;

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common
