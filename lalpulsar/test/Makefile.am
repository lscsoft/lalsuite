TEST_CPPFLAGS = -DTEST_DATA_DIR='"$(abs_srcdir)/"'
LDADD = $(top_builddir)/src/liblalpulsar.la

check_PROGRAMS = \
	LALBarycenterTest \
	BinarySSBTimesTest \
	ExtrapolatePulsarSpinsTest \
	FitToPulsarTest \
	FlatMeshTest \
	GeneralMeshTest \
	GeneralMetricTest \
	GeneratePulsarSignalTest \
	GetAMCoeffsTest \
	LatticeTilingTest \
	LineRobustStatsTest \
	NewGetAMCoeffsTest \
	XLALComputeAMTest \
	PtoleMeshTest \
	PtoleMetricTest \
	SFTfileIOTest \
	SimulateTaylorCWTest \
	SkyMetricTest \
	StackMetricTest \
	TestConstructPLUT \
	TestDriveHough \
	TestDriveNDHough \
	TestHoughMap \
	LFTandTSutilsTest \
	TestNDConstructPLUT \
	TestNDHoughMap \
	TestNDPeak2PHMD \
	TestPeak2PHMD \
	TestStatistics \
	TestVelocity \
	TwoDMeshTest \
	UniversalDopplerMetricTest \
	NormalizeSFTRngMedTest \
	XLALMultiNoiseWeightsTest \
	FileIOTest \
	ComputeFstatTest
#	SFTCleanTest

if LALXML
check_PROGRAMS += LALPulsarXMLTest
endif

TESTS = $(check_PROGRAMS)
EXTRA_PROGRAMS = XLALLoadSFTsTest

TESTS_ENVIRONMENT = \
	LAL_DEBUG_LEVEL=msglvl1,memdbg; export LAL_DEBUG_LEVEL;

TwoDMeshTest_SOURCES = TwoDMeshTest.c TwoDMeshPlot.h TwoDMeshPlot.c

UniversalDopplerMetricTest_SOURCES = UniversalDopplerMetricTest.c old-FstatMetric.c

noinst_DATA = \
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
	MultiNoiseWeightsTest1.sft \
	MultiNoiseWeightsTest2.sft \
	MultiNoiseWeightsTest3.sft \
	MultiNoiseWeightsTest4.sft \
	MultiNoiseWeightsTest5.sft \
	MultiNoiseWeightsTest6.sft \
	MultiNoiseWeightsTest7.sft \
	inputsft.0 \
	inputsft.1 \
	earth98.dat \
	sun98.dat \
	testTimestamps.dat

MOSTLYCLEANFILES = \
	H-*_H1*.sft \
	OutHistogram.asc \
	OutHough.asc \
	outputsft*.sft \
	TS_R4.dat \
	LFT_R4.dat \
	LFT_C8.dat

pkgdata_DATA = \
	circularEphem.dat \
	earth00-19-DE200.dat.gz \
	earth00-19-DE405.dat.gz \
	earth00-19-DE414.dat.gz \
	sun00-19-DE200.dat.gz \
	sun00-19-DE405.dat.gz \
	sun00-19-DE414.dat.gz \
	testTimestamps.dat \
	VOTable-1.1.xsd \
	tdb_2000-2019.dat.gz \
	te405_2000-2019.dat.gz

EXTRA_DIST = $(DATA)
