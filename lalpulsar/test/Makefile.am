AM_CFLAGS = $(LALSUPPORT_CFLAGS)
INCLUDES = -I../include

if LALSUPPORT
check_PROGRAMS = \
	ExtrapolatePulsarSpinsTest \
	FitToPulsarTest \
	FoldAmplitudesTest \
	GeneralMeshTest \
	GeneralMetricTest \
	GeneratePulsarSignalTest \
	GetAMCoeffsTest \
	NewGetAMCoeffsTest \
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
	UniversalDopplerMetricTest

##	SFTCleanTest \
##	NormalizeSFTRngMedTest \
else
check_PROGRAMS =
endif

ExtrapolatePulsarSpinsTest_SOURCES = ExtrapolatePulsarSpinsTest.c
ExtrapolatePulsarSpinsTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

FitToPulsarTest_SOURCES = FitToPulsarTest.c
FitToPulsarTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

FoldAmplitudesTest_SOURCES = FoldAmplitudesTest.c
FoldAmplitudesTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

GeneralMeshTest_SOURCES = GeneralMeshTest.c
GeneralMeshTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

GeneralMetricTest_SOURCES = GeneralMetricTest.c
GeneralMetricTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

GeneratePulsarSignalTest_SOURCES = GeneratePulsarSignalTest.c
GeneratePulsarSignalTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

GetAMCoeffsTest_SOURCES = GetAMCoeffsTest.c
GetAMCoeffsTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

NewGetAMCoeffsTest_SOURCES = NewGetAMCoeffsTest.c
NewGetAMCoeffsTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

#NormalizeSFTRngMedTest_SOURCES = NormalizeSFTRngMedTest.c
#NormalizeSFTRngMedTest_LDADD = \
#	$(LALSUPPORT_LIBS) \
#	$(top_builddir)/src/liblalpulsar.la

PtoleMeshTest_SOURCES = PtoleMeshTest.c
PtoleMeshTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

PtoleMetricTest_SOURCES = PtoleMetricTest.c
PtoleMetricTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

PulsarCatTest_SOURCES = PulsarCatTest.c
PulsarCatTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

ResampleTest_SOURCES = ResampleTest.c
ResampleTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

#SFTCleanTest_SOURCES = SFTCleanTest.c
#SFTCleanTest_LDADD = \
#	$(LALSUPPORT_LIBS) \
#	$(top_builddir)/src/liblalpulsar.la

SFTfileIOTest_SOURCES = SFTfileIOTest.c
SFTfileIOTest_CPPFLAGS = -DSRCDIR=\"$(srcdir)\"
SFTfileIOTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

SimulateTaylorCWTest_SOURCES = SimulateTaylorCWTest.c
SimulateTaylorCWTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

SkyMetricTest_SOURCES = SkyMetricTest.c
SkyMetricTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

StackMetricTest_SOURCES = StackMetricTest.c
StackMetricTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

TestConstructPLUT_SOURCES = TestConstructPLUT.c
TestConstructPLUT_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

TestDriveHough_SOURCES = TestDriveHough.c
TestDriveHough_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

TestDriveNDHough_SOURCES = TestDriveNDHough.c
TestDriveNDHough_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

TestHoughMap_SOURCES = TestHoughMap.c
TestHoughMap_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

TestNDConstructPLUT_SOURCES = TestNDConstructPLUT.c
TestNDConstructPLUT_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

TestNDHoughMap_SOURCES = TestNDHoughMap.c
TestNDHoughMap_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

TestNDPeak2PHMD_SOURCES = TestNDPeak2PHMD.c
TestNDPeak2PHMD_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

TestPeak2PHMD_SOURCES = TestPeak2PHMD.c
TestPeak2PHMD_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

TestStatistics_SOURCES = TestStatistics.c
TestStatistics_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

TestVelocity_SOURCES = TestVelocity.c
TestVelocity_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

TwoDMeshTest_SOURCES = TwoDMeshTest.c TwoDMeshPlot.h TwoDMeshPlot.c
TwoDMeshTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

UniversalDopplerMetricTest_SOURCES = UniversalDopplerMetricTest.c
UniversalDopplerMetricTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalpulsar.la

if LALSUPPORT
TESTS = $(check_PROGRAMS)
else
TESTS =
endif

TESTS_ENVIRONMENT = 1>/dev/null LAL_DATA_PATH=$(srcdir)

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
	inputsft.1

MOSTLYCLEANFILES = \
	H-1_H1*.sft \
	OutHistogram.asc \
	OutHough.asc \
	outputsft*.sft

pkgdata_DATA = \
	circularEphem.dat \
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
	earth98.dat \
	earth99.dat \
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
	sun98.dat \
	sun99.dat


EXTRA_DIST = $(DATA)
