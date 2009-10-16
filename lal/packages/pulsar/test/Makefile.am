## Process this file with automake to produce Makefile.in

## global settings for all programs built here
LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

check_PROGRAMS = StackMetricTest FlatMeshTest TwoDMeshTest SkyMetricTest DirectedMeshTest ResampleTest \
	FoldAmplitudesTest LALBarycenterTest PtoleMetricTest PtoleMeshTest GeneralMetricTest GeneralMeshTest \
	FitToPulsarTest PulsarCatTest ExtrapolatePulsarSpinsTest GetAMCoeffsTest NewGetAMCoeffsTest

StackMetricTest_SOURCES = StackMetricTest.c
FlatMeshTest_SOURCES = FlatMeshTest.c
TwoDMeshTest_SOURCES = TwoDMeshTest.c TwoDMeshPlot.h TwoDMeshPlot.c
SkyMetricTest_SOURCES = SkyMetricTest.c
DirectedMeshTest_SOURCES = DirectedMeshTest.c
ResampleTest_SOURCES = ResampleTest.c
##LALDemodTest_SOURCES = LALDemodTest.c
FoldAmplitudesTest_SOURCES = FoldAmplitudesTest.c
LALBarycenterTest_SOURCES = LALBarycenterTest.c
PtoleMetricTest_SOURCES = PtoleMetricTest.c
PtoleMeshTest_SOURCES = PtoleMeshTest.c
GeneralMetricTest_SOURCES = GeneralMetricTest.c
GeneralMeshTest_SOURCES = GeneralMeshTest.c
FitToPulsarTest_SOURCES = FitToPulsarTest.c
PulsarCatTest_SOURCES = PulsarCatTest.c
ExtrapolatePulsarSpinsTest_SOURCES = ExtrapolatePulsarSpinsTest.c
GetAMCoeffsTest_SOURCES = GetAMCoeffsTest.c
NewGetAMCoeffsTest_SOURCES = NewGetAMCoeffsTest.c


pkgdata_DATA = eigen_Tf0-1e8_n-2.dat eigenInv_Tf0-1e8_n-2.dat \
  range_tauOverT-1e5_n-2.dat in.data \
  earth98.dat sun98.dat \
  earth99.dat sun99.dat \
  earth00.dat sun00.dat \
  earth01.dat sun01.dat \
  earth02.dat sun02.dat \
  earth03.dat sun03.dat \
  earth04.dat sun04.dat \
  earth05.dat sun05.dat  \
  earth06.dat sun06.dat  \
  earth07.dat sun07.dat  \
  earth08.dat sun08.dat  \
  earth09.dat sun09.dat  \
  earth00-04.dat sun00-04.dat \
  earth03-06.dat sun03-06.dat \
  earth05-09.dat sun05-09.dat \
  earth09-11.dat sun09-11.dat \
  earth09-11-DE405.dat sun09-11-DE405.dat \
  earth98.ilwd sun98.ilwd \
  earth99.ilwd sun99.ilwd \
  earth00.ilwd sun00.ilwd \
  earth01.ilwd sun01.ilwd \
  earth02.ilwd sun02.ilwd \
  earth03.ilwd sun03.ilwd \
  earth04.ilwd sun04.ilwd \
  earth00-04.ilwd sun00-04.ilwd

EXTRA_DIST = $(DATA) ephemtoilwd.sh
EXTRA_TEXSOURCES =
MOSTLYCLEANFILES = .dvi-dep times.blah


TESTS = StackMetricTest FlatMeshTest TwoDMeshTest SkyMetricTest      \
  ResampleTest FoldAmplitudesTest LALBarycenterTest \
  PtoleMetricTest PtoleMeshTest GeneralMetricTest GeneralMeshTest FitToPulsarTest PulsarCatTest \
  ExtrapolatePulsarSpinsTest GetAMCoeffsTest NewGetAMCoeffsTest


TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null LAL_DATA_PATH=$(srcdir)
DISTCLEANFILES = *.tex
MAINTAINERCLEANFILES = earth*.ilwd sun*.ilwd
dvi-local: .dvi-dep

SUFFIXES=.ilwd .dat
.dat.ilwd:; $(srcdir)/ephemtoilwd.sh $< > $(srcdir)/$@
@DVIMKFILE@
