BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

SFTReferenceLibrary = SFTReferenceLibrary.c SFTReferenceLibrary.h

bin_PROGRAMS = \
	lalapps_ComputePSD \
	lalapps_GetSFTVelocities \
	lalapps_SFTclean \
	lalapps_SFTvalidate  \
	lalapps_SFTwrite \
	lalapps_compareSFTs \
	lalapps_dumpSFT \
	lalapps_splitSFTs \
	$(END_OF_LIST)

lalapps_ComputePSD_SOURCES = ComputePSD.c
lalapps_GetSFTVelocities_SOURCES = GetSFTVelocities.c
lalapps_SFTclean_SOURCES = SFTclean.c
lalapps_SFTvalidate_SOURCES = SFTvalidate.c $(SFTReferenceLibrary)
lalapps_SFTwrite_SOURCES = SFTwrite.c $(SFTReferenceLibrary)
lalapps_compareSFTs_SOURCES = compareSFTs.c
lalapps_dumpSFT_SOURCES = dumpSFT.c
lalapps_splitSFTs_SOURCES = splitSFTs.c $(SFTReferenceLibrary)

EXTRA_DIST = \
	SFTSpecification.dox \
	testComputePSD.sh \
	testSFTvalidate.sh \
	$(END_OF_LIST)

CLEANFILES = \
	SFT-bad1[0-4] \
	SFT-bad[1-9] \
	SFT-good \
	SFT-test[1-7] \
	$(END_OF_LIST)

TESTS = \
	testComputePSD.sh \
	testSFTvalidate.sh \
	$(END_OF_LIST)
