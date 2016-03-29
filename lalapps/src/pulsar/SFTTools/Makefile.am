BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/lalapps.am

SFTReferenceLibrary = SFTReferenceLibrary.c SFTReferenceLibrary.h

bin_PROGRAMS = \
	lalapps_ConvertToSFTv2 \
	lalapps_SFTdumpall \
	lalapps_SFTdumpheader \
	lalapps_SFTvalidate  \
	lalapps_SFTwrite \
	lalapps_splitSFTs \
	$(END_OF_LIST)

lalapps_ConvertToSFTv2_SOURCES = ConvertToSFTv2.c
lalapps_SFTdumpall_SOURCES = SFTdumpall.c $(SFTReferenceLibrary)
lalapps_SFTdumpheader_SOURCES = SFTdumpheader.c $(SFTReferenceLibrary)
lalapps_SFTvalidate_SOURCES = SFTvalidate.c $(SFTReferenceLibrary)
lalapps_SFTwrite_SOURCES = SFTwrite.c $(SFTReferenceLibrary)
lalapps_splitSFTs_SOURCES = splitSFTs.c $(SFTReferenceLibrary)

EXTRA_DIST = \
	SFTReferenceLibrary.dox \
	testSFTvalidate.sh \
	$(END_OF_LIST)

CLEANFILES = \
	SFT-bad1[0-4] \
	SFT-bad[1-9] \
	SFT-good \
	SFT-test[1-7] \
	$(END_OF_LIST)

TESTS_ENVIRONMENT = \
	LAL_DEBUG_LEVEL=msglvl1,memdbg; export LAL_DEBUG_LEVEL; \
	LAL_DATA_PATH="$(LAL_DATA_PATH)"; export LAL_DATA_PATH;

TESTS = \
	testSFTvalidate.sh \
	$(END_OF_LIST)
