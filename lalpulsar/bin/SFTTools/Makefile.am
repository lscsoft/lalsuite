BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps_pulsar.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

bin_PROGRAMS = \
	lalpulsar_ComputePSD \
	lalpulsar_SFTclean \
	lalpulsar_SFTvalidate  \
	lalpulsar_WriteSFTsfromSFDBs \
	lalpulsar_compareSFTs \
	lalpulsar_dumpSFT \
	lalpulsar_splitSFTs \
	$(END_OF_LIST)

lalpulsar_ComputePSD_SOURCES = \
	ComputePSD.c \
	$(END_OF_LIST)

lalpulsar_SFTclean_SOURCES = \
	SFTclean.c \
	$(END_OF_LIST)

lalpulsar_SFTvalidate_SOURCES = \
	SFTvalidate.c \
	$(END_OF_LIST)

lalpulsar_compareSFTs_SOURCES = \
	compareSFTs.c \
	$(END_OF_LIST)

lalpulsar_dumpSFT_SOURCES = \
	dumpSFT.c \
	$(END_OF_LIST)

lalpulsar_splitSFTs_SOURCES = \
	splitSFTs.c \
	$(END_OF_LIST)

lalpulsar_WriteSFTsfromSFDBs_SOURCES = \
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
