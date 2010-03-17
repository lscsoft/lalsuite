AM_CFLAGS = $(LALSUPPORT_CFLAGS)
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalinspiral.la

if LALSUPPORT
check_PROGRAMS = \
	BasicInjectTest \
	CoarseTest \
	FindChirpBCVSpinTest \
	FindChirpBankVetoTest \
	FindChirpTDTest \
	GeneratePPNAmpCorInspiralTest \
	GeneratePPNInspiralTest \
	InjectionInterfaceTest \
	InspiralBCVSpinBankTest \
	InspiralSpinBankTest \
	LALInspiralSpinningBHBinariesTest \
	LALInspiralTest \
	LALSTPNWaveformTest \
	MetricTest \
	RandomInspiralSignalTest
TESTS = $(check_PROGRAMS)
else
check_PROGRAMS =
TESTS =
endif

# non-building test
# FilterTest
# FindChirpACTDTest
# FindChirpBCVTest
# FindChirpSPTest
# FindChirpTmpltTest
# InspiralSpinBankwNDTemplateBankTest
# NRWaveInjectTest
# RingTest

# failing test
# NRWaveIOTest

# requires lalstochastic
# CoherentInspiralFilterTest
# CoherentInspiralInputTest

# requires lalframe
# LALHybridTest

EXTRA_DIST = \
	NRWaveIOTest.cfg \
	NRWaveIOTest.data \
	RandomInspiralSignalTest.par \
	injection.xml

CLEANFILES = \
	*.out \
	*.dat

TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null
