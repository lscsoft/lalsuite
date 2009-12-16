LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null LAL_DATA_PATH=$(srcdir)

TESTS = \
	ExtrapolatePulsarSpinsTest \
	FitToPulsarTest \
	FlatMeshTest \
	FoldAmplitudesTest \
	GeneralMeshTest \
	GeneralMetricTest \
	GetAMCoeffsTest \
	NewGetAMCoeffsTest \
	PtoleMeshTest \
	PtoleMetricTest \
	PulsarCatTest \
	ResampleTest \
	SkyMetricTest \
	StackMetricTest \
	TwoDMeshTest \
	UniversalDopplerMetricTest

check_PROGRAMS = \
	DirectedMeshTest \
	ExtrapolatePulsarSpinsTest \
	FitToPulsarTest \
	FlatMeshTest \
	FoldAmplitudesTest \
	GeneralMeshTest \
	GeneralMetricTest \
	GetAMCoeffsTest \
	NewGetAMCoeffsTest \
	PtoleMeshTest \
	PtoleMetricTest \
	PulsarCatTest \
	ResampleTest \
	SkyMetricTest \
	StackMetricTest \
	TwoDMeshTest \
	UniversalDopplerMetricTest

DirectedMeshTest_SOURCES = DirectedMeshTest.c
ExtrapolatePulsarSpinsTest_SOURCES = ExtrapolatePulsarSpinsTest.c
FitToPulsarTest_SOURCES = FitToPulsarTest.c
FlatMeshTest_SOURCES = FlatMeshTest.c
FoldAmplitudesTest_SOURCES = FoldAmplitudesTest.c
GeneralMeshTest_SOURCES = GeneralMeshTest.c
GeneralMetricTest_SOURCES = GeneralMetricTest.c
GetAMCoeffsTest_SOURCES = GetAMCoeffsTest.c
#LALDemodTest_SOURCES = LALDemodTest.c
NewGetAMCoeffsTest_SOURCES = NewGetAMCoeffsTest.c
PtoleMeshTest_SOURCES = PtoleMeshTest.c
PtoleMetricTest_SOURCES = PtoleMetricTest.c
PulsarCatTest_SOURCES = PulsarCatTest.c
ResampleTest_SOURCES = ResampleTest.c
SkyMetricTest_SOURCES = SkyMetricTest.c
StackMetricTest_SOURCES = StackMetricTest.c
TwoDMeshTest_SOURCES = TwoDMeshTest.c TwoDMeshPlot.h TwoDMeshPlot.c
UniversalDopplerMetricTest_SOURCES = UniversalDopplerMetricTest.c

pkgdata_DATA = \
	circularEphem.dat \
	earth00-04.dat \
	earth00-04.ilwd \
	earth00.dat \
	earth00.ilwd \
	earth01.dat \
	earth01.ilwd \
	earth02.dat \
	earth02.ilwd \
	earth03-06.dat \
	earth03.dat \
	earth03.ilwd \
	earth04.dat \
	earth04.ilwd \
	earth05-09.dat \
	earth05.dat \
	earth06.dat \
	earth07.dat \
	earth08.dat \
	earth09-11-DE405.dat \
	earth09-11.dat \
	earth09.dat \
	earth98.dat \
	earth98.ilwd \
	earth99.dat \
	earth99.ilwd \
	eigenInv_Tf0-1e8_n-2.dat \
	eigen_Tf0-1e8_n-2.dat \
	in.data \
	range_tauOverT-1e5_n-2.dat \
	sun00-04.dat \
	sun00-04.ilwd \
	sun00.dat \
	sun00.ilwd \
	sun01.dat \
	sun01.ilwd \
	sun02.dat \
	sun02.ilwd \
	sun03-06.dat \
	sun03.dat \
	sun03.ilwd \
	sun04.dat \
	sun04.ilwd \
	sun05-09.dat \
	sun05.dat  \
	sun06.dat  \
	sun07.dat  \
	sun08.dat  \
	sun09-11-DE405.dat \
	sun09-11.dat \
	sun09.dat  \
	sun98.dat \
	sun98.ilwd \
	sun99.dat \
	sun99.ilwd

EXTRA_TEXSOURCES =

EXTRA_DIST = $(DATA) ephemtoilwd.sh

MOSTLYCLEANFILES = .dvi-dep times.blah
DISTCLEANFILES = *.tex
MAINTAINERCLEANFILES = earth*.ilwd sun*.ilwd

dvi-local: .dvi-dep

SUFFIXES=.ilwd .dat
.dat.ilwd:; $(srcdir)/ephemtoilwd.sh $< > $(srcdir)/$@

@DVIMKFILE@
