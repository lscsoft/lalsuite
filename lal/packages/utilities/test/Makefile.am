## Process this file with automake to produce Makefile.in
check_PROGRAMS = RandomTest MersenneRandomTest FindRootTest \
  IntegrateTest InterpolateTest SortTest ODETest DirichletTest \
  DetInverseTest EigenTest LALRunningMedianTest RngMedBiasTest \
  SCoarseGrainFrequencySeriesTest DCoarseGrainFrequencySeriesTest \
  CCoarseGrainFrequencySeriesTest ZCoarseGrainFrequencySeriesTest \
  SFTCleanTest
RandomTest_SOURCES = RandomTest.c
RandomTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
MersenneRandomTest_SOURCES = MersenneRandomTest.c
MersenneRandomTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
FindRootTest_SOURCES = FindRootTest.c
FindRootTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
IntegrateTest_SOURCES = IntegrateTest.c
IntegrateTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
InterpolateTest_SOURCES = InterpolateTest.c
InterpolateTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
SortTest_SOURCES = SortTest.c
SortTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
ODETest_SOURCES = ODETest.c
ODETest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
DirichletTest_SOURCES = DirichletTest.c
DirichletTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
SCoarseGrainFrequencySeriesTest_SOURCES = SCoarseGrainFrequencySeriesTest.c CheckStatus.c
SCoarseGrainFrequencySeriesTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
DCoarseGrainFrequencySeriesTest_SOURCES = DCoarseGrainFrequencySeriesTest.c CheckStatus.c
DCoarseGrainFrequencySeriesTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
CCoarseGrainFrequencySeriesTest_SOURCES = CCoarseGrainFrequencySeriesTest.c CheckStatus.c
CCoarseGrainFrequencySeriesTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
ZCoarseGrainFrequencySeriesTest_SOURCES = ZCoarseGrainFrequencySeriesTest.c CheckStatus.c
ZCoarseGrainFrequencySeriesTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
DetInverseTest_SOURCES = DetInverseTest.c
DetInverseTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
EigenTest_SOURCES = EigenTest.c
EigenTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALRunningMedianTest_SOURCES = LALRunningMedianTest.c
LALRunningMedianTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
RngMedBiasTest_SOURCES = RngMedBiasTest.c
RngMedBiasTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
SFTCleanTest_SOURCES = SFTCleanTest.c
SFTCleanTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
EXTRA_DIST = CheckStatus.c CheckStatus.h
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = RandomTest MersenneRandomTest FindRootTest IntegrateTest \
  InterpolateTest SortTest ODETest DirichletTest EigenTest \
  DetInverseTest LALRunningMedianTest RngMedBiasTest \
  SCoarseGrainFrequencySeriesTest  DCoarseGrainFrequencySeriesTest \
  CCoarseGrainFrequencySeriesTest  ZCoarseGrainFrequencySeriesTest
TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null
MOSTLYCLEANFILES = .dvi-dep PrintVector.* *.out
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
