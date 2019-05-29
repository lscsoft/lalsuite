EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_test.am

# Add compiled test programs to this variable
test_programs += BinarySSBTimesTest
test_programs += ComputeFstatTest
test_programs += ConstructPLUTTest
test_programs += CWSignalBandTest
test_programs += DopplerScanTest
test_programs += DriveHoughTest
test_programs += DriveNDHoughTest
test_programs += ExtrapolatePulsarSpinsTest
test_programs += FITSFileIOTest
test_programs += FileIOTest
test_programs += GeneralMeshTest
test_programs += GeneralMetricTest
test_programs += GeneratePulsarSignalTest
test_programs += HoughMapTest
test_programs += LALBarycenterTest
test_programs += LFTandTSutilsTest
test_programs += LatticeTilingTest
test_programs += LineRobustStatsTest
test_programs += MetricUtilsTest
test_programs += NDConstructPLUTTest
test_programs += NDHoughMapTest
test_programs += NDPeak2PHMDTest
test_programs += NormalizeSFTRngMedTest
test_programs += Peak2PHMDTest
test_programs += PtoleMeshTest
test_programs += PtoleMetricTest
test_programs += ReadTEMPOFileTest
test_programs += SFTfileIOTest
test_programs += SimulateTaylorCWTest
test_programs += StatisticsTest
test_programs += SuperskyMetricsTest
test_programs += TwoDMeshTest
test_programs += UniversalDopplerMetricTest
test_programs += VelocityTest
test_programs += XLALComputeAMTest
test_programs += XLALMultiNoiseWeightsTest

# Add shell, Python, etc. test scripts to this variable
test_scripts += pulsar_toa_test.sh
if SWIG_BUILD_PYTHON
test_scripts += simulateCWTest.py
test_scripts += heterodyned_model_test.py
endif

# Add any helper programs required by tests to this variable
test_helpers += TEMPOcomparison
test_helpers += PulsarTOATest

TwoDMeshTest_SOURCES = \
	TwoDMeshPlot.c \
	TwoDMeshPlot.h \
	TwoDMeshTest.c \
	$(END_OF_LIST)

UniversalDopplerMetricTest_SOURCES = \
	UniversalDopplerMetricTest.c \
	old-FstatMetric.c \
	$(END_OF_LIST)

MOSTLYCLEANFILES = \
	FITSFileIOTest.fits \
	H-*_H1*.sft \
	LFT_C8.dat \
	LFT_R4.dat \
	LatticeTilingTest.fits \
	OutHistogram.asc \
	OutHough.asc \
	SuperskyMetricsTest.fits \
	TEMPOcomparison.par \
	TEMPOcomparison.tim \
	TS_R4.dat \
	V-1_V1_1800SFT_simCW_simulateCWTest-*.sft \
	outputsft*.sft \
	$(END_OF_LIST)

EXTRA_DIST += \
	MultiNoiseWeightsTest1.sft \
	MultiNoiseWeightsTest2.sft \
	MultiNoiseWeightsTest3.sft \
	MultiNoiseWeightsTest4.sft \
	MultiNoiseWeightsTest5.sft \
	MultiNoiseWeightsTest6.sft \
	MultiNoiseWeightsTest7.sft \
	SFT-bad1 \
	SFT-bad10 \
	SFT-bad11 \
	SFT-bad12 \
	SFT-bad13 \
	SFT-bad14 \
	SFT-bad2 \
	SFT-bad3 \
	SFT-bad4 \
	SFT-bad5 \
	SFT-bad6 \
	SFT-bad7 \
	SFT-bad8 \
	SFT-bad9 \
	SFT-good \
	SFT-test1 \
	SFT-test2 \
	SFT-test3 \
	SFT-test4 \
	SFT-test5 \
	SFT-test6 \
	SFT-test7 \
	circularEphem.dat \
	earth98.dat \
	inputsft.0 \
	inputsft.1 \
	sun98.dat \
	testTimestamps.dat \
	testTimestampsNew.dat \
	TOAtest.par \
	TOAtest.simulate \
	pks2gps.clk \
	$(END_OF_LIST)
