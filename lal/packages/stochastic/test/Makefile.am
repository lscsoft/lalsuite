# Process this file with automake to produce Makefile.in
check_PROGRAMS =  StochasticCrossCorrelationTest  SZeroPadAndFFTTest \
                  StochasticOptimalFilterTest     StochasticOmegaGWTest \
                  OverlapReductionFunctionTest    StochasticInverseNoiseTest \
                  DirichletTest

#                   SFalseAlarmTest       SFalseDismissalTest  \
#                   DFalseAlarmTest       SFDismissalVsFAlarmTest \
#                   DFalseDismissalTest   DFDismissalVsFAlarmTest

StochasticCrossCorrelationTest_SOURCES = StochasticCrossCorrelationTest.c \
                                         CheckStatus.c
StochasticCrossCorrelationTest_LDADD = $(top_builddir)/lib/liblal.la
SZeroPadAndFFTTest_SOURCES = SZeroPadAndFFTTest.c CheckStatus.c
SZeroPadAndFFTTest_LDADD = $(top_builddir)/lib/liblal.la
StochasticOptimalFilterTest_SOURCES = StochasticOptimalFilterTest.c \
                                         CheckStatus.c
StochasticOptimalFilterTest_LDADD = $(top_builddir)/lib/liblal.la
StochasticOmegaGWTest_SOURCES = StochasticOmegaGWTest.c CheckStatus.c
StochasticOmegaGWTest_LDADD = $(top_builddir)/lib/liblal.la
OverlapReductionFunctionTest_SOURCES = OverlapReductionFunctionTest.c \
                                         CheckStatus.c
OverlapReductionFunctionTest_LDADD = $(top_builddir)/lib/liblal.la
StochasticInverseNoiseTest_SOURCES = StochasticInverseNoiseTest.c \
                                         CheckStatus.c
StochasticInverseNoiseTest_LDADD = $(top_builddir)/lib/liblal.la
# SFalseAlarmTest_SOURCES = SFalseAlarmTest.c
# SFalseAlarmTest_LDADD = $(top_builddir)/lib/liblal.la
# SFalseDismissalTest_SOURCES = SFalseDismissalTest.c
# SFalseDismissalTest_LDADD = $(top_builddir)/lib/liblal.la
# SFDismissalVsFAlarmTest_SOURCES = SFDismissalVsFAlarmTest.c
# SFDismissalVsFAlarmTest_LDADD = $(top_builddir)/lib/liblal.la
# DFalseAlarmTest_SOURCES = DFalseAlarmTest.c
# DFalseAlarmTest_LDADD = $(top_builddir)/lib/liblal.la
# DFalseDismissalTest_SOURCES = DFalseDismissalTest.c
# DFalseDismissalTest_LDADD = $(top_builddir)/lib/liblal.la
# DFDismissalVsFAlarmTest_SOURCES = DFDismissalVsFAlarmTest.c
# DFDismissalVsFAlarmTest_LDADD = $(top_builddir)/lib/liblal.la
DirichletTest_SOURCES = DirichletTest.c
DirichletTest_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

TESTS =	  StochasticCrossCorrelationTest         SZeroPadAndFFTTest \
          StochasticOptimalFilterTest           StochasticOmegaGWTest \
          OverlapReductionFunctionTest           StochasticInverseNoiseTest \
          DirichletTest

#    SFalseAlarmTest       SFalseDismissalTest  \
#           DFalseAlarmTest       SFDismissalVsFAlarmTest \
#           DFalseDismissalTest   DFDismissalVsFAlarmTest \

MOSTLYCLEANFILES = *.dat* PrintVector.0*

EXTRA_DIST = CheckStatus.c CheckStatus.h
