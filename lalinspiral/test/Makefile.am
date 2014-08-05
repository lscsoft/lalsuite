EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_test.am

# Add compiled test programs to this variable
test_programs += BCVTemplates
test_programs += BasicInjectTest
test_programs += CoarseTest
test_programs += CoarseTest2
test_programs += FindChirpBankVetoTest
test_programs += GenerateInspiralWaveform
test_programs += GeneratePPNAmpCorInspiralTest
test_programs += GeneratePPNInspiralTest
test_programs += GetOrientationEllipse
test_programs += InjectionInterfaceTest
test_programs += InspiralBCVSpinBankTest
test_programs += InspiralSpinBankTest
test_programs += LALInspiralSpinningBHBinariesTest
test_programs += LALInspiralTaylorT2Test
test_programs += LALInspiralTaylorT3Test
test_programs += LALInspiralTaylorT4Test
test_programs += LALInspiralTest
test_programs += LALSTPNWaveformTest
test_programs += MetricTest
test_programs += MetricTestBCV
test_programs += MetricTestPTF
test_programs += PNTemplates
test_programs += RandomInspiralSignalTest
# non-building tests:
#test_programs += BCVSpinTemplates
#test_programs += ChirpSpace
#test_programs += FilterTest
#test_programs += FindChirpACTDTest
#test_programs += FindChirpBCVTest
#test_programs += FindChirpSPTest
#test_programs += FindChirpTmpltTest
#test_programs += InspiralSpinBankwNDTemplateBankTest
#test_programs += NRWaveInjectTest
#test_programs += RingTest
#test_programs += TestRingdown
# failing tests:
#test_programs += FindChirpTDTest
#test_programs += LALHybridTest
#test_programs += NRWaveIOTest
#test_programs += SpaceCovering
# requires lalstochastic:
#test_programs += CoherentInspiralFilterTest
#test_programs += CoherentInspiralInputTest
# requires lalpulsar:
#test_programs += BCVTemplatesFlatMesh

# Add shell, Python, etc. test scripts to this variable
test_scripts +=

# Add any helper programs required by tests to this variable
test_helpers +=

MOSTLYCLEANFILES = \
	*.dat \
	*.out \
	$(END_OF_LIST)

EXTRA_DIST += \
	NRWaveIOTest.cfg \
	NRWaveIOTest.data \
	RandomInspiralSignalTest.par \
	injection.xml \
	$(END_OF_LIST)
