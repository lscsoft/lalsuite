include $(top_srcdir)/gnuscripts/lalsuite_header_links.am

TEST_CPPFLAGS = -DTEST_DATA_DIR='"$(abs_srcdir)/"'

check_PROGRAMS = \
	CCoarseGrainFrequencySeriesTest \
	CSInterpolateTest \
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

CSInterpolateTest_SOURCES = \
	CSInterpolateTest.c

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

MOSTLYCLEANFILES = \
	PrintVector.* \
	*.out
