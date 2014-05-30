TEST_CPPFLAGS = -DTEST_DATA_DIR='"$(abs_srcdir)/"'
LDADD = $(top_builddir)/src/liblalstochastic.la

check_PROGRAMS = \
	StochasticCrossCorrelationStatisticTest \
	StochasticHeterodynedCrossCorrelationStatisticTest \
	StochasticCrossCorrelationSpectrumTest \
	SZeroPadAndFFTTest \
	CZeroPadAndFFTTest \
	StochasticOptimalFilterTest \
	StochasticOmegaGWTest \
	StochasticOptimalFilterNormalizationTest \
	OverlapReductionFunctionTest \
	StochasticInverseNoiseTest \
	AstroOmegaTest \
	SimulateSBTest \
	SimulatePopcornTest
TESTS = $(check_PROGRAMS)

StochasticCrossCorrelationStatisticTest_SOURCES = \
  StochasticCrossCorrelationStatisticTest.c CheckStatus.c
StochasticHeterodynedCrossCorrelationStatisticTest_SOURCES = \
  StochasticHeterodynedCrossCorrelationStatisticTest.c CheckStatus.c
StochasticCrossCorrelationSpectrumTest_SOURCES = \
  StochasticCrossCorrelationSpectrumTest.c CheckStatus.c
SZeroPadAndFFTTest_SOURCES = SZeroPadAndFFTTest.c CheckStatus.c
CZeroPadAndFFTTest_SOURCES = CZeroPadAndFFTTest.c CheckStatus.c
StochasticOptimalFilterTest_SOURCES = \
  StochasticOptimalFilterTest.c CheckStatus.c
StochasticOptimalFilterNormalizationTest_SOURCES = \
  StochasticOptimalFilterNormalizationTest.c CheckStatus.c
StochasticOmegaGWTest_SOURCES = \
  StochasticOmegaGWTest.c CheckStatus.c
OverlapReductionFunctionTest_SOURCES = \
  OverlapReductionFunctionTest.c CheckStatus.c
StochasticInverseNoiseTest_SOURCES = \
  StochasticInverseNoiseTest.c CheckStatus.c
SimulateSBTest_SOURCES = SimulateSBTest.c CheckStatus.c

EXTRA_DIST = \
	pipeline \
	spectrum1.in \
	spectrum2.in \
	response1.in \
	response2.in \
	CheckStatus.c \
	CheckStatus.h \
	L1.dat \
	H2.dat

MOSTLYCLEANFILES = \
	??_*.dat \
	PrintVector.0* \
	WHITENED-SB*
