AM_CFLAGS = $(LALSUPPORT_CFLAGS) -DDATADIR='"$(abs_srcdir)/"'
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalpulsar.la

if LALXML
XMLTEST = LALPulsarXMLTest
endif

check_PROGRAMS = \
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
	ResampleTest \
	SFTfileIOTest \
	SimulateTaylorCWTest \
	SkyMetricTest \
	StackMetricTest \
	TestConstructPLUT \
	TestDriveHough \
	TestDriveNDHough \
	TestHoughMap \
	TestNDConstructPLUT \
	TestNDHoughMap \
	TestNDPeak2PHMD \
	TestPeak2PHMD \
	TestStatistics \
	TestVelocity \
	TwoDMeshTest \
	UniversalDopplerMetricTest \
	NormalizeSFTRngMedTest \
	$(XMLTEST)
##	SFTCleanTest \

TESTS = $(check_PROGRAMS)
EXTRA_PROGRAMS = XLALLoadSFTsTest

AM_CPPFLAGS = -DSRCDIR=\"$(srcdir)\"
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
	inputsft.0 \
	inputsft.1 \
	testTimestamps.dat

MOSTLYCLEANFILES = \
	H-1_H1*.sft \
	OutHistogram.asc \
	OutHough.asc \
	outputsft*.sft

pkgdata_DATA = \
	circularEphem.dat \
	earth00-19-DE200.dat.gz \
	earth00-19-DE405.dat.gz \
	earth00-19-DE414.dat.gz \
	earth00-04.dat \
	earth00.dat \
	earth01.dat \
	earth02.dat \
	earth03-06.dat \
	earth03.dat \
	earth04.dat \
	earth05-09.dat \
	earth05.dat \
	earth06.dat \
	earth07.dat \
	earth08.dat \
	earth09-11-DE405.dat \
	earth09-11.dat \
	earth09.dat \
	earth10.dat \
	earth11.dat \
	earth12.dat \
	earth13.dat \
	earth10-13.dat \
	earth98.dat \
	earth99.dat \
	sun00-19-DE200.dat.gz \
        sun00-19-DE405.dat.gz \
        sun00-19-DE414.dat.gz \
	sun00-04.dat \
	sun00.dat \
	sun01.dat \
	sun02.dat \
	sun03-06.dat \
	sun03.dat \
	sun04.dat \
	sun05-09.dat \
	sun05.dat  \
	sun06.dat  \
	sun07.dat  \
	sun08.dat  \
	sun09-11-DE405.dat \
	sun09-11.dat \
	sun09.dat  \
	sun10.dat \
	sun11.dat \
	sun12.dat \
	sun13.dat \
	sun10-13.dat \
	sun98.dat \
	sun99.dat \
	testTimestamps.dat \
	VOTable-1.1.xsd \
	tdb_2000-2019.dat.gz \
	te405_2000-2019.dat.gz

EXTRA_DIST = $(DATA)
