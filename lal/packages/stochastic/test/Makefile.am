# Process this file with automake to produce Makefile.in
check_PROGRAMS = 	  StochasticCrossCorrelationStatisticTest \
	  StochasticHeterodynedCrossCorrelationStatisticTest \
          StochasticCrossCorrelationSpectrumTest \
	  SZeroPadAndFFTTest                    CZeroPadAndFFTTest \
          StochasticOptimalFilterTest           StochasticOmegaGWTest \
          StochasticOptimalFilterNormalizationTest \
          OverlapReductionFunctionTest          StochasticInverseNoiseTest \
          DirichletTest        \
	  SCoarseGrainFrequencySeriesTest    DCoarseGrainFrequencySeriesTest \
	  CCoarseGrainFrequencySeriesTest    ZCoarseGrainFrequencySeriesTest

#                   SFalseAlarmTest       SFalseDismissalTest  \
#                   DFalseAlarmTest       SFDismissalVsFAlarmTest \
#                   DFalseDismissalTest   DFDismissalVsFAlarmTest

StochasticCrossCorrelationStatisticTest_SOURCES = \
                               StochasticCrossCorrelationStatisticTest.c \
                               CheckStatus.c
StochasticCrossCorrelationStatisticTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
StochasticHeterodynedCrossCorrelationStatisticTest_SOURCES = \
                    StochasticHeterodynedCrossCorrelationStatisticTest.c \
                               CheckStatus.c
StochasticHeterodynedCrossCorrelationStatisticTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
StochasticCrossCorrelationSpectrumTest_SOURCES = \
                   StochasticCrossCorrelationSpectrumTest.c \
                               CheckStatus.c
StochasticCrossCorrelationSpectrumTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
SZeroPadAndFFTTest_SOURCES = SZeroPadAndFFTTest.c CheckStatus.c
SZeroPadAndFFTTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
CZeroPadAndFFTTest_SOURCES = CZeroPadAndFFTTest.c CheckStatus.c
CZeroPadAndFFTTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
StochasticOptimalFilterTest_SOURCES = StochasticOptimalFilterTest.c \
                                         CheckStatus.c
StochasticOptimalFilterTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
StochasticOptimalFilterNormalizationTest_SOURCES = StochasticOptimalFilterNormalizationTest.c \
                                         CheckStatus.c
StochasticOptimalFilterNormalizationTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
StochasticOmegaGWTest_SOURCES = StochasticOmegaGWTest.c CheckStatus.c
StochasticOmegaGWTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
OverlapReductionFunctionTest_SOURCES = OverlapReductionFunctionTest.c \
                                         CheckStatus.c
OverlapReductionFunctionTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
StochasticInverseNoiseTest_SOURCES = StochasticInverseNoiseTest.c \
                                         CheckStatus.c
StochasticInverseNoiseTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
# SFalseAlarmTest_SOURCES = SFalseAlarmTest.c
# SFalseAlarmTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
# SFalseDismissalTest_SOURCES = SFalseDismissalTest.c
# SFalseDismissalTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
# SFDismissalVsFAlarmTest_SOURCES = SFDismissalVsFAlarmTest.c
# SFDismissalVsFAlarmTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
# DFalseAlarmTest_SOURCES = DFalseAlarmTest.c
# DFalseAlarmTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
# DFalseDismissalTest_SOURCES = DFalseDismissalTest.c
# DFalseDismissalTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
# DFDismissalVsFAlarmTest_SOURCES = DFDismissalVsFAlarmTest.c
# DFDismissalVsFAlarmTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
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
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

TESTS =	  StochasticCrossCorrelationStatisticTest \
	  StochasticHeterodynedCrossCorrelationStatisticTest \
          StochasticCrossCorrelationSpectrumTest \
	  SZeroPadAndFFTTest                    CZeroPadAndFFTTest \
          StochasticOptimalFilterTest           \
          StochasticOptimalFilterNormalizationTest   StochasticOmegaGWTest \
          OverlapReductionFunctionTest           StochasticInverseNoiseTest \
          DirichletTest        \
	  SCoarseGrainFrequencySeriesTest  DCoarseGrainFrequencySeriesTest \
	  CCoarseGrainFrequencySeriesTest  ZCoarseGrainFrequencySeriesTest

#    SFalseAlarmTest       SFalseDismissalTest  \
#           DFalseAlarmTest       SFDismissalVsFAlarmTest \
#           DFalseDismissalTest   DFDismissalVsFAlarmTest

EXTRA_DIST = CheckStatus.c CheckStatus.h pipeline spectrum1.in spectrum2.in \
             response1.in response2.in

MOSTLYCLEANFILES = .dvi-dep *.dat* PrintVector.0*
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
