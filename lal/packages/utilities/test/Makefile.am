AM_CPPFLAGS = -I$(top_builddir)/include -I$(srcdir)/../include @AM_CPPFLAGS@
LDADD = \
	$(top_builddir)/lib/lal/liblal.la \
	$(top_builddir)/lib/lalsupport/src/liblalsupport.la

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
