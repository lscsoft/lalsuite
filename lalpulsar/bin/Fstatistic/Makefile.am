BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps_pulsar.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

bin_PROGRAMS = \
	lalpulsar_ComputeFstatBenchmark \
	lalpulsar_ComputeFstatLatticeCount \
	lalpulsar_ComputeFstatMCUpperLimit \
	lalpulsar_ComputeFstatistic_v2 \
	lalpulsar_PredictFstat \
	lalpulsar_compareFstats\
	lalpulsar_synthesizeBstatMC \
	lalpulsar_synthesizeLVStats \
	lalpulsar_synthesizeTransientStats \
	$(END_OF_LIST)

lalpulsar_compareFstats_SOURCES = \
	compareFstats.c \
	$(END_OF_LIST)

lalpulsar_ComputeFstatistic_v2_SOURCES = \
	ComputeFstatistic_v2.c \
	$(END_OF_LIST)

lalpulsar_PredictFstat_SOURCES = \
	PredictFstat.c \
	$(END_OF_LIST)

lalpulsar_ComputeFstatMCUpperLimit_SOURCES = \
	ComputeFstatMCUpperLimit.c \
	$(END_OF_LIST)

lalpulsar_ComputeFstatLatticeCount_SOURCES = \
	ComputeFstatLatticeCount.c \
	$(END_OF_LIST)

lalpulsar_synthesizeBstatMC_SOURCES = \
	synthesizeBstatMC.c \
	$(END_OF_LIST)

lalpulsar_synthesizeTransientStats_SOURCES = \
	synthesizeTransientStats.c \
	$(END_OF_LIST)

lalpulsar_synthesizeLVStats_SOURCES = \
	synthesizeLVStats.c \
	$(END_OF_LIST)

lalpulsar_ComputeFstatBenchmark_SOURCES = \
	ComputeFstatBenchmark.c \
	$(END_OF_LIST)

# Add shell test scripts to this variable
test_scripts += testPredictFstat.sh
test_scripts += testComputeFstatistic_v2.sh
test_scripts += testComputeFstatistic_v2_singleF.sh
test_scripts += testComputeFstatistic_v2_binary.sh
test_scripts += testComputeFstatistic_v2_grids.sh
test_scripts += testComputeFstatistic_v2_resamp.sh
test_scripts += testComputeFstatistic_v2_transient.sh
test_scripts += testComputeFstatBenchmark.sh
test_scripts += testComputeFstatMCUpperLimit.sh
test_scripts += test_synthesizeBstatMC.sh
test_scripts += test_synthesizeLVStats.sh
test_scripts += test_synthesizeTransientStats.sh

# Add any helper programs required by tests to this variable
test_helpers += SemiAnalyticF

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=
