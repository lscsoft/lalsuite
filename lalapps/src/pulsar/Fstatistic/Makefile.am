MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

bin_PROGRAMS = \
	lalapps_ComputeFstatBenchmark \
	lalapps_ComputeFstatLatticeCount \
	lalapps_ComputeFstatMCUpperLimit \
	lalapps_ComputeFstatistic_v2 \
	lalapps_PredictFstat \
	lalapps_SemiAnalyticF \
	lalapps_compareFstats\
	lalapps_synthesizeBstatMC \
	lalapps_synthesizeLVStats \
	lalapps_synthesizeTransientStats \
	$(END_OF_LIST)

lalapps_SemiAnalyticF_SOURCES = SemiAnalyticF.c

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

EXTRA_DIST = \
	testCFSv2.sh \
	testCFSv2_binary.sh \
	testCFSv2_grids.sh \
	testCFSv2_resamp.sh \
	testCFSv2_transient.sh \
	testPredictFstat.sh \
	testSynthLV.sh \
	$(END_OF_LIST)

EXTRA_DIST += \
	H1_test_timestamps.dat \
	L1_test_timestamps.dat \
	testCFSv2.dat.ref.gz \
	testCFSv2_grid0.dat.ref.gz \
	testCFSv2_grid1.dat.ref.gz \
	testCFSv2_grid2.dat.ref.gz \
	testCFSv2_grid3.dat.ref.gz \
	testCFSv2_grid6.dat.ref.gz \
	testCFSv2_grid8.dat.ref.gz \
	testCFSv2_grid9.dat.ref.gz \
	$(END_OF_LIST)

TESTS = \
	testPredictFstat.sh \
	testSynthLV.sh \
	$(END_OF_LIST)

if !CONDOR_ENABLED
TESTS += \
	testCFSv2.sh \
	testCFSv2_binary.sh \
	testCFSv2_grids.sh \
	testCFSv2_resamp.sh \
	testCFSv2_transient.sh \
	$(END_OF_LIST)
endif
