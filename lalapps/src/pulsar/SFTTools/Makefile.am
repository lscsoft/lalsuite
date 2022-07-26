BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps_pulsar.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

bin_PROGRAMS = \
	lalapps_ComputePSD \
	lalapps_SFTclean \
	lalapps_SFTvalidate  \
	lalapps_WriteSFTsfromSFDBs \
	lalapps_compareSFTs \
	lalapps_dumpSFT \
	lalapps_splitSFTs \
	$(END_OF_LIST)

lalapps_ComputePSD_SOURCES = \
	ComputePSD.c \
	$(END_OF_LIST)

lalapps_SFTclean_SOURCES = \
	SFTclean.c \
	$(END_OF_LIST)

lalapps_SFTvalidate_SOURCES = \
	SFTvalidate.c \
	$(END_OF_LIST)

lalapps_compareSFTs_SOURCES = \
	compareSFTs.c \
	$(END_OF_LIST)

lalapps_dumpSFT_SOURCES = \
	dumpSFT.c \
	$(END_OF_LIST)

lalapps_splitSFTs_SOURCES = \
	splitSFTs.c \
	$(END_OF_LIST)

lalapps_WriteSFTsfromSFDBs_SOURCES = \
	WriteSFTsfromSFDBs.c \
	$(END_OF_LIST)

# Add shell test scripts to this variable
test_scripts += testComputePSD.sh
test_scripts += testSFTvalidate.sh
test_scripts += testdumpSFT.sh
test_scripts += testcompareSFTs.sh
test_scripts += testsplitSFTs.sh
test_scripts += testSFTclean.sh
if HAVE_PYTHON
test_scripts += testWriteSFTsfromSFDBs.py
endif

# Add any helper programs required by tests to this variable
test_helpers += SFTwrite

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=
