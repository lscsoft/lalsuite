## Process this file with automake to produce Makefile.in
check_PROGRAMS = StackMetricTest FlatMeshTest StreamInputTest ResampleTest LALDemodTest
noinst_DATA = eigen_Tf0-1e8_n-2.dat eigenInv_Tf0-1e8_n-2.dat range_tauOverT-1e5_n-2.dat StreamInput.dat in.data
EXTRA_DIST = $(DATA)
MOSTLYCLEANFILES = .dvi-dep times.blah
StackMetricTest_SOURCES = StackMetricTest.c
StackMetricTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
FlatMeshTest_SOURCES = FlatMeshTest.c StreamVectorInput.c StreamVectorSequenceInput.c
FlatMeshTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
StreamInputTest_SOURCES = StreamInputTest.c StreamInput.h StreamVectorInput.c StreamVectorSequenceInput.c
StreamInputTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
ResampleTest_SOURCES = ResampleTest.c
ResampleTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
LALDemodTest_SOURCES = LALDemodTest.c
LALDemodTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = StackMetricTest FlatMeshTest StreamInputTest ResampleTest
TESTS_ENVIRONMENT = LAL_DATA_PATH=$(srcdir)
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
