## Process this file with automake to produce Makefile.in
check_PROGRAMS = StackMetricTest FlatMeshTest TwoDMeshTest SkyMetricTest DirectedMeshTest ResampleTest FoldAmplitudesTest LALBarycenterTest PtoleMetricTest PtoleMeshTest GeneralMetricTest GeneralMeshTest HeterodynePulsarTest FitToPulsarTest PulsarCatTest
#check_PROGRAMS = StackMetricTest FlatMeshTest TwoDMeshTest ResampleTest LALDemodTest FoldAmplitudesTest LALBarycenterTest PtoleMetricTest PtoleMeshTest GeneralMetricTest GeneralMeshTest HeterodynePulsarTest FitToPulsarTest PulsarCatTest
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
  earth98.ilwd sun98.ilwd \
  earth99.ilwd sun99.ilwd \
  earth00.ilwd sun00.ilwd \
  earth01.ilwd sun01.ilwd \
  earth02.ilwd sun02.ilwd \
  earth03.ilwd sun03.ilwd \
  earth04.ilwd sun04.ilwd \
  earth00-04.ilwd sun00-04.ilwd
#EXTRA_DIST = $(DATA) ephemtoilwd.sh
EXTRA_DIST = $(DATA) LALDemodTest.c ephemtoilwd.sh
EXTRA_TEXSOURCES = LALDemodTest.c
MOSTLYCLEANFILES = .dvi-dep times.blah
StackMetricTest_SOURCES = StackMetricTest.c
StackMetricTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
FlatMeshTest_SOURCES = FlatMeshTest.c
FlatMeshTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TwoDMeshTest_SOURCES = TwoDMeshTest.c TwoDMeshPlot.h TwoDMeshPlot.c
TwoDMeshTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
SkyMetricTest_SOURCES = SkyMetricTest.c
SkyMetricTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
DirectedMeshTest_SOURCES = DirectedMeshTest.c
DirectedMeshTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
ResampleTest_SOURCES = ResampleTest.c
ResampleTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
#LALDemodTest_SOURCES = LALDemodTest.c
#LALDemodTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
FoldAmplitudesTest_SOURCES = FoldAmplitudesTest.c
FoldAmplitudesTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALBarycenterTest_SOURCES = LALBarycenterTest.c
LALBarycenterTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
PtoleMetricTest_SOURCES = PtoleMetricTest.c
PtoleMetricTest_LDADD = $(top_builddir)/lib/liblal.la
PtoleMeshTest_SOURCES = PtoleMeshTest.c
PtoleMeshTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
GeneralMetricTest_SOURCES = GeneralMetricTest.c
GeneralMetricTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
GeneralMeshTest_SOURCES = GeneralMeshTest.c
GeneralMeshTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
HeterodynePulsarTest_SOURCES = HeterodynePulsarTest.c
HeterodynePulsarTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
FitToPulsarTest_SOURCES = FitToPulsarTest.c
FitToPulsarTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
PulsarCatTest_SOURCES = PulsarCatTest.c
PulsarCatTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
#TESTS = StackMetricTest FlatMeshTest TwoDMeshTest ResampleTest FoldAmplitudesTest LALBarycenterTest PtoleMetricTest PtoleMeshTest GeneralMetricTest GeneralMeshTest FitToPulsarTest LALDemodTest
TESTS = StackMetricTest FlatMeshTest TwoDMeshTest SkyMetricTest      \
  ResampleTest FoldAmplitudesTest LALBarycenterTest \
  PtoleMetricTest PtoleMeshTest GeneralMetricTest GeneralMeshTest FitToPulsarTest PulsarCatTest
TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null LAL_DATA_PATH=$(srcdir)
DISTCLEANFILES = *.tex
MAINTAINERCLEANFILES = earth*.ilwd sun*.ilwd
dvi-local: .dvi-dep

SUFFIXES=.ilwd .dat
.dat.ilwd:; $(srcdir)/ephemtoilwd.sh $< > $(srcdir)/$@
@DVIMKFILE@
