BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalapps_pulsar_test.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

bin_PROGRAMS = \
	lalapps_ComputeFstatBenchmark \
	lalapps_ComputeFstatLatticeCount \
	lalapps_ComputeFstatMCUpperLimit \
	lalapps_ComputeFstatistic_v2 \
	lalapps_PredictFstat \
	lalapps_compareFstats\
	lalapps_synthesizeBstatMC \
	lalapps_synthesizeLVStats \
	lalapps_synthesizeTransientStats \
	$(END_OF_LIST)

lalapps_compareFstats_SOURCES = compareFstats.c

lalapps_ComputeFstatistic_v2_SOURCES = \
	ComputeFstatistic_v2.c \
	HeapToplist.c \
	HeapToplist.h \
	$(END_OF_LIST)

lalapps_PredictFstat_SOURCES = PredictFstat.c

lalapps_ComputeFstatMCUpperLimit_SOURCES = ComputeFstatMCUpperLimit.c

lalapps_ComputeFstatLatticeCount_SOURCES = ComputeFstatLatticeCount.c

lalapps_synthesizeBstatMC_SOURCES = synthesizeBstatMC.c

lalapps_synthesizeTransientStats_SOURCES = synthesizeTransientStats.c

lalapps_synthesizeLVStats_SOURCES = synthesizeLVStats.c

lalapps_ComputeFstatBenchmark_SOURCES = ComputeFstatBenchmark.c

# Add shell test scripts to this variable
test_scripts += testPredictFstat.sh
test_scripts += testComputeFstatistic_v2.sh
test_scripts += testComputeFstatistic_v2_singleF.sh
test_scripts += testComputeFstatistic_v2_binary.sh
test_scripts += testComputeFstatistic_v2_grids.sh
test_scripts += testComputeFstatistic_v2_resamp.sh
test_scripts += testComputeFstatistic_v2_transient.sh
test_scripts += test_synthesizeLVStats.sh

# Add any helper programs required by tests to this variable
test_helpers += SemiAnalyticF

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=
