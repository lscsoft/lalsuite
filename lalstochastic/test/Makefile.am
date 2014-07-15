EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_test.am

# Add compiled test programs to this variable
test_programs += AstroOmegaTest
test_programs += CZeroPadAndFFTTest
test_programs += OverlapReductionFunctionTest
test_programs += SZeroPadAndFFTTest
test_programs += SimulatePopcornTest
test_programs += SimulateSBTest
test_programs += StochasticCrossCorrelationSpectrumTest
test_programs += StochasticCrossCorrelationStatisticTest
test_programs += StochasticHeterodynedCrossCorrelationStatisticTest
test_programs += StochasticInverseNoiseTest
test_programs += StochasticOmegaGWTest
test_programs += StochasticOptimalFilterNormalizationTest
test_programs += StochasticOptimalFilterTest

# Add shell, Python, etc. test scripts to this variable
test_scripts +=

# Add any helper programs required by tests to this variable
test_helpers +=

MOSTLYCLEANFILES = \
	??_*.dat \
	PrintVector.0* \
	WHITENED-SB* \
	$(END_OF_LIST)

EXTRA_DIST += \
	CheckStatus.h \
	H2.dat \
	L1.dat \
	pipeline \
	response1.in \
	response2.in \
	spectrum1.in \
	spectrum2.in \
	$(END_OF_LIST)
