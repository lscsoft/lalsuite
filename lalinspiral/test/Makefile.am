AM_CFLAGS = $(LALSUPPORT_CFLAGS)
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalinspiral.la

if LALSUPPORT
check_PROGRAMS = \
	BCVTemplates \
	BCVTemplatesFlatMesh \
	BasicInjectTest \
	CoarseTest \
	CoarseTest2 \
	FindChirpBCVSpinTest \
	FindChirpBankVetoTest \
	FindChirpTDTest \
	GenerateInspiralWaveform \
	GeneratePPNAmpCorInspiralTest \
	GeneratePPNInspiralTest \
	GetOrientationEllipse \
	InjectionInterfaceTest \
	InspiralBCVSpinBankTest \
	InspiralSpinBankTest \
	LALInspiralSpinningBHBinariesTest \
	LALInspiralTest \
	LALSTPNWaveformTest \
	MetricTest \
	MetricTestBCV \
	MetricTestPTF \
	PNTemplates \
	RandomInspiralSignalTest \
	SpaceCovering
TESTS = $(check_PROGRAMS)
endif

# non-building test
# BCVSpinTemplates
# ChirpSpace
# FilterTest
# FindChirpACTDTest
# FindChirpBCVTest
# FindChirpSPTest
# FindChirpTmpltTest
# InspiralSpinBankwNDTemplateBankTest
# NRWaveInjectTest
# RingTest
# TestRingdown

# failing test
# LALHybridTest
# NRWaveIOTest

# requires lalstochastic
# CoherentInspiralFilterTest
# CoherentInspiralInputTest

EXTRA_DIST = \
	NRWaveIOTest.cfg \
	NRWaveIOTest.data \
	RandomInspiralSignalTest.par \
	injection.xml

CLEANFILES = \
	*.out \
	*.dat

TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null
