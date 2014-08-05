EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_test.am

# Add compiled test programs to this variable
test_programs += BinarySSBTimesTest
test_programs += ComputeFstatTest
test_programs += ConstructPLUTTest
test_programs += DriveHoughTest
test_programs += DriveNDHoughTest
test_programs += ExtrapolatePulsarSpinsTest
test_programs += FileIOTest
test_programs += FitToPulsarTest
test_programs += FlatMeshTest
test_programs += GeneralMeshTest
test_programs += GeneralMetricTest
test_programs += GeneratePulsarSignalTest
test_programs += GetAMCoeffsTest
test_programs += HoughMapTest
test_programs += LALBarycenterTest
test_programs += LFTandTSutilsTest
test_programs += LatticeTilingTest
test_programs += LineRobustStatsTest
test_programs += NDConstructPLUTTest
test_programs += NDHoughMapTest
test_programs += NDPeak2PHMDTest
test_programs += NewGetAMCoeffsTest
test_programs += NormalizeSFTRngMedTest
test_programs += Peak2PHMDTest
test_programs += PtoleMeshTest
test_programs += PtoleMetricTest
test_programs += SFTfileIOTest
test_programs += SimulateTaylorCWTest
test_programs += SkyMetricTest
test_programs += StackMetricTest
test_programs += StatisticsTest
test_programs += TwoDMeshTest
test_programs += UniversalDopplerMetricTest
test_programs += VelocityTest
test_programs += XLALComputeAMTest
test_programs += XLALMultiNoiseWeightsTest
#test_programs += XLALLoadSFTsTest
#test_programs += SFTCleanTest
if LALXML
test_programs += LALPulsarXMLTest
endif

# Add shell, Python, etc. test scripts to this variable
test_scripts +=

# Add any helper programs required by tests to this variable
test_helpers +=

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
	H-*_H1*.sft \
	LFT_C8.dat \
	LFT_R4.dat \
	OutHistogram.asc \
	OutHough.asc \
	TS_R4.dat \
	outputsft*.sft \
	$(END_OF_LIST)

pkgdata_DATA = \
	VOTable-1.1.xsd \
	circularEphem.dat \
	earth00-19-DE200.dat.gz \
	earth00-19-DE405.dat.gz \
	earth00-19-DE414.dat.gz \
	sun00-19-DE200.dat.gz \
	sun00-19-DE405.dat.gz \
	sun00-19-DE414.dat.gz \
	tdb_2000-2019.dat.gz \
	te405_2000-2019.dat.gz \
	testTimestamps.dat \
	$(END_OF_LIST)

EXTRA_DIST += \
	$(pkgdata_DATA) \
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
	earth98.dat \
	inputsft.0 \
	inputsft.1 \
	sun98.dat \
	testTimestamps.dat \
	$(END_OF_LIST)
