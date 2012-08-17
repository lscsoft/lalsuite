INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

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
	CheckStatus.h \
	CCoarseGrainFrequencySeriesTest.c

DCoarseGrainFrequencySeriesTest_SOURCES = \
	CheckStatus.c \
	CheckStatus.h \
	DCoarseGrainFrequencySeriesTest.c

SCoarseGrainFrequencySeriesTest_SOURCES = \
	CheckStatus.c \
	CheckStatus.h \
	SCoarseGrainFrequencySeriesTest.c

ZCoarseGrainFrequencySeriesTest_SOURCES = \
	CheckStatus.c \
	CheckStatus.h \
	ZCoarseGrainFrequencySeriesTest.c

EXTRA_DIST = \
	CheckStatus.c \
	CheckStatus.h

TESTS = \
	$(check_PROGRAMS)

TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null

MOSTLYCLEANFILES = \
	PrintVector.* \
	*.out
