# Process this file with automake to produce Makefile.in

AMCFLAGS = $(LALSUPPORT_CFLAGS)
LIBS := $(LIBS) $(LALSUPPORT_LIBS)

if LALSUPPORT
check_PROGRAMS = StochasticCrossCorrelationStatisticTest \
  StochasticHeterodynedCrossCorrelationStatisticTest \
  StochasticCrossCorrelationSpectrumTest SZeroPadAndFFTTest \
  CZeroPadAndFFTTest StochasticOptimalFilterTest StochasticOmegaGWTest \
  StochasticOptimalFilterNormalizationTest \
  OverlapReductionFunctionTest StochasticInverseNoiseTest \
  AstroOmegaTest SimulateSBTest SimulatePopcornTest
else
check_PROGRAMS =
endif

StochasticCrossCorrelationStatisticTest_SOURCES = \
  StochasticCrossCorrelationStatisticTest.c CheckStatus.c
StochasticCrossCorrelationStatisticTest_LDADD = \
  $(top_builddir)/src/liblalstochastic.la
StochasticHeterodynedCrossCorrelationStatisticTest_SOURCES = \
  StochasticHeterodynedCrossCorrelationStatisticTest.c CheckStatus.c
StochasticHeterodynedCrossCorrelationStatisticTest_LDADD = \
  $(top_builddir)/src/liblalstochastic.la
StochasticCrossCorrelationSpectrumTest_SOURCES = \
  StochasticCrossCorrelationSpectrumTest.c CheckStatus.c
StochasticCrossCorrelationSpectrumTest_LDADD = \
  $(top_builddir)/src/liblalstochastic.la
SZeroPadAndFFTTest_SOURCES = SZeroPadAndFFTTest.c CheckStatus.c
SZeroPadAndFFTTest_LDADD = \
  $(top_builddir)/src/liblalstochastic.la
CZeroPadAndFFTTest_SOURCES = CZeroPadAndFFTTest.c CheckStatus.c
CZeroPadAndFFTTest_LDADD = \
  $(top_builddir)/src/liblalstochastic.la
StochasticOptimalFilterTest_SOURCES = \
  StochasticOptimalFilterTest.c CheckStatus.c
StochasticOptimalFilterTest_LDADD = \
  $(top_builddir)/src/liblalstochastic.la
StochasticOptimalFilterNormalizationTest_SOURCES = \
  StochasticOptimalFilterNormalizationTest.c CheckStatus.c
StochasticOptimalFilterNormalizationTest_LDADD = \
  $(top_builddir)/src/liblalstochastic.la
StochasticOmegaGWTest_SOURCES = \
  StochasticOmegaGWTest.c CheckStatus.c
StochasticOmegaGWTest_LDADD = \
  $(top_builddir)/src/liblalstochastic.la
OverlapReductionFunctionTest_SOURCES = \
  OverlapReductionFunctionTest.c CheckStatus.c
OverlapReductionFunctionTest_LDADD = \
  $(top_builddir)/src/liblalstochastic.la
StochasticInverseNoiseTest_SOURCES = \
  StochasticInverseNoiseTest.c CheckStatus.c
StochasticInverseNoiseTest_LDADD = \
  $(top_builddir)/src/liblalstochastic.la
SimulateSBTest_SOURCES = SimulateSBTest.c CheckStatus.c
SimulateSBTest_LDADD = \
  $(top_builddir)/src/liblalstochastic.la
SimulatePopcornTest_SOURCES = SimulatePopcornTest.c
SimulatePopcornTest_LDADD = \
  $(top_builddir)/src/liblalstochastic.la
AstroOmegaTest_SOURCES = AstroOmegaTest.c
AstroOmegaTest_LDADD = \
  $(top_builddir)/src/liblalstochastic.la

INCLUDES = -I../include

if LALSUPPORT
TESTS = $(check_PROGRAMS)
else
TESTS =
endif

TESTS_ENVIRONMENT = 1>/dev/null

EXTRA_DIST = pipeline spectrum1.in spectrum2.in response1.in \
  response2.in CheckStatus.c CheckStatus.h L1.dat H2.dat

MOSTLYCLEANFILES = ??_*.dat PrintVector.0* WHITENED-SB*
