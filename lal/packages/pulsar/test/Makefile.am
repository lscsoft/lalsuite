## Process this file with automake to produce Makefile.in
#check_PROGRAMS = StackMetricTest FlatMeshTest TwoDMeshTest ResampleTest LALDemodTest FoldAmplitudesTest LALBarycenterTest PtoleMetricTest
check_PROGRAMS = StackMetricTest FlatMeshTest TwoDMeshTest ResampleTest FoldAmplitudesTest LALBarycenterTest PtoleMetricTest
noinst_DATA = eigen_Tf0-1e8_n-2.dat eigenInv_Tf0-1e8_n-2.dat range_tauOverT-1e5_n-2.dat in.data earth98.dat sun98.dat
#EXTRA_DIST = $(DATA)
EXTRA_DIST = $(DATA) LALDemodTest.c
EXTRA_TEXSOURCES = LALDemodTest.c
MOSTLYCLEANFILES = .dvi-dep times.blah
StackMetricTest_SOURCES = StackMetricTest.c
StackMetricTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
FlatMeshTest_SOURCES = FlatMeshTest.c
FlatMeshTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TwoDMeshTest_SOURCES = TwoDMeshTest.c TwoDMeshPlot.h TwoDMeshPlot.c
TwoDMeshTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
ResampleTest_SOURCES = ResampleTest.c
ResampleTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
#LALDemodTest_SOURCES = LALDemodTest.c
#LALDemodTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
FoldAmplitudesTest_SOURCES = FoldAmplitudesTest.c
FoldAmplitudesTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALBarycenterTest_SOURCES = LALBarycenterTest.c LALInitBarycenter.h LALInitBarycenter.c
LALBarycenterTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
PtoleMetricTest_SOURCES = PtoleMetricTest.c
PtoleMetricTest_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = StackMetricTest FlatMeshTest TwoDMeshTest ResampleTest FoldAmplitudesTest LALBarycenterTest PtoleMetricTest
TESTS_ENVIRONMENT = LAL_DATA_PATH=$(srcdir)
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
