AM_CFLAGS = @AM_CFLAGS@ $(LALSUPPORT_CFLAGS) -DTEST_DATA_DIR='"$(abs_srcdir)/"'
AM_CPPFLAGS = @AM_CPPFLAGS@ -I$(top_builddir)/include -DSRCDIR=\"$(srcdir)\"
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalpulsar.la

if LALXML
XMLTEST = LALPulsarXMLTest
endif

check_PROGRAMS = \
	LALBarycenterTest \
	BinarySSBTimesTest \
	ExtrapolatePulsarSpinsTest \
	FitToPulsarTest \
	FlatMeshTest \
	FoldAmplitudesTest \
	GeneralMeshTest \
	GeneralMetricTest \
	GeneratePulsarSignalTest \
	GetAMCoeffsTest \
	NewGetAMCoeffsTest \
	XLALComputeAMTest \
	PtoleMeshTest \
	PtoleMetricTest \
	PulsarCatTest \
	SFTfileIOTest \
	SimulateTaylorCWTest \
	SkyMetricTest \
	StackMetricTest \
	TestConstructPLUT \
	TestDriveHough \
	TestDriveNDHough \
	TestHoughMap \
	TestLFTandTSutils \
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
	$(XMLTEST)
##	SFTCleanTest \

TESTS = $(check_PROGRAMS)
EXTRA_PROGRAMS = XLALLoadSFTsTest

TESTS_ENVIRONMENT = LAL_DEBUG_LEVEL=memdbg

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
	testLFTandTSutils-LFT.sft \
	testLFTandTSutils-timeseries.dat

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
