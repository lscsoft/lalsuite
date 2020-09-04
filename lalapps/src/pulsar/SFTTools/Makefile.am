BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalapps_pulsar_test.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

bin_PROGRAMS = \
	lalapps_ComputePSD \
	lalapps_SFTclean \
	lalapps_SFTvalidate  \
	lalapps_compareSFTs \
	lalapps_dumpSFT \
	lalapps_splitSFTs \
	$(END_OF_LIST)

lalapps_ComputePSD_SOURCES = ComputePSD.c
lalapps_SFTclean_SOURCES = SFTclean.c
lalapps_SFTvalidate_SOURCES = SFTvalidate.c
lalapps_compareSFTs_SOURCES = compareSFTs.c
lalapps_dumpSFT_SOURCES = dumpSFT.c
lalapps_splitSFTs_SOURCES = splitSFTs.c

EXTRA_DIST += \
	SFTSpecification.dox \
	$(END_OF_LIST)

# Add shell test scripts to this variable
test_scripts += testComputePSD.sh
test_scripts += testSFTvalidate.sh
test_scripts += testdumpSFT.sh
test_scripts += testcompareSFTs.sh
test_scripts += testsplitSFTs.sh

# Add any helper programs required by tests to this variable
test_helpers += SFTwrite

SFTwrite_SOURCES = SFTwrite.c
