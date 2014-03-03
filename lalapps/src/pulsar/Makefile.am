## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(top_builddir)/src/lalapps -DPKG_DATA_DIR='"$(pkgdatadir)/"' @AM_CPPFLAGS@
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

SUBDIRS = \
	MakeSFTs \
	Injections \
	FDS_isolated \
	hough \
	GCT \
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
	lalapps_FreqAverager \
	lalapps_psd \
	lalapps_ComputePSD \
	lalapps_FreqAverager_v2 \
	lalapps_GetSFTVelocities \
	lalapps_PrintDetectorState \
	lalapps_ComputeAntennaPattern

lalapps_FreqAverager_SOURCES = FreqAverager.c
lalapps_psd_SOURCES = psd.c psd.h
lalapps_ComputePSD_SOURCES = ComputePSD.c
lalapps_FreqAverager_v2_SOURCES =  FreqAverager_v2.c
lalapps_GetSFTVelocities_SOURCES = GetSFTVelocities.c
lalapps_PrintDetectorState_SOURCES = PrintDetectorState.c
lalapps_ComputeAntennaPattern_SOURCES = ComputeAntennaPattern.c

EXTRA_DIST = \
	FreqAverager.c \
	psd.c \
	psd.h \
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
