INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

check_PROGRAMS = \
	CCoarseGrainFrequencySeriesTest \
	DCoarseGrainFrequencySeriesTest \
	DetInverseTest \
	DirichletTest \
	EigenTest \
	FindRootTest \
	IntegrateTest \
	InterpolateTest \
	LALRunningMedianTest \
	MersenneRandomTest \
	ODETest \
	RandomTest \
	RngMedBiasTest \
	SCoarseGrainFrequencySeriesTest \
	SortTest \
	ZCoarseGrainFrequencySeriesTest

CCoarseGrainFrequencySeriesTest_SOURCES = \
	CheckStatus.c \
	CCoarseGrainFrequencySeriesTest.c
CCoarseGrainFrequencySeriesTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

DCoarseGrainFrequencySeriesTest_SOURCES = \
	CheckStatus.c \
	DCoarseGrainFrequencySeriesTest.c
DCoarseGrainFrequencySeriesTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

DetInverseTest_SOURCES = DetInverseTest.c
DetInverseTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

DirichletTest_SOURCES = DirichletTest.c
DirichletTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

EigenTest_SOURCES = EigenTest.c
EigenTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

FindRootTest_SOURCES = FindRootTest.c
FindRootTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

IntegrateTest_SOURCES = IntegrateTest.c
IntegrateTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

InterpolateTest_SOURCES = InterpolateTest.c
InterpolateTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

LALRunningMedianTest_SOURCES = LALRunningMedianTest.c
LALRunningMedianTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

MersenneRandomTest_SOURCES = MersenneRandomTest.c
MersenneRandomTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

ODETest_SOURCES = ODETest.c
ODETest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

RandomTest_SOURCES = RandomTest.c
RandomTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

SCoarseGrainFrequencySeriesTest_SOURCES = \
	CheckStatus.c \
	SCoarseGrainFrequencySeriesTest.c
SCoarseGrainFrequencySeriesTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

SortTest_SOURCES = SortTest.c
SortTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

RngMedBiasTest_SOURCES = RngMedBiasTest.c
RngMedBiasTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

ZCoarseGrainFrequencySeriesTest_SOURCES = \
	CheckStatus.c \
	ZCoarseGrainFrequencySeriesTest.c
ZCoarseGrainFrequencySeriesTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

EXTRA_DIST = \
	CheckStatus.c \
	CheckStatus.h

TESTS = \
	CCoarseGrainFrequencySeriesTest \
	DCoarseGrainFrequencySeriesTest \
	DetInverseTest \
	DirichletTest \
	EigenTest \
	FindRootTest \
	IntegrateTest \
	InterpolateTest \
	LALRunningMedianTest \
	MersenneRandomTest \
	ODETest \
	RandomTest \
	RngMedBiasTest \
	SCoarseGrainFrequencySeriesTest \
	SortTest \
	ZCoarseGrainFrequencySeriesTest

TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null

MOSTLYCLEANFILES = .dvi-dep PrintVector.* *.out
DISTCLEANFILES = *.tex

dvi-local: .dvi-dep

@DVIMKFILE@
